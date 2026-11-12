"""Reconstruct Drassodes integrated PM-eye log contrast; see docs/spider_prc_reconstruction.md.

This is a reconstruction from recovered source and investigator instructions,
not a certified reproduction of historical site observations.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
import os
import platform
from pathlib import Path
import subprocess

import numpy as np


VERSION = "spider-eye-reconstruction-2"
ROOT = Path(__file__).resolve().parents[1]
GROUPS = [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def new_directory(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=False)
    return path


def relative_path(value, base):
    p = Path(value).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def discover(roots, output, max_dirs=100000):
    """Inventory without opening image pixels, following symlinks or changing inputs."""
    candidates, other_image_folders, archives, warnings, seen = [], [], [], [], set()
    count = 0
    for root in roots:
        root = Path(root).expanduser().resolve()
        if not root.is_dir():
            warnings.append(f"Missing search root: {root}")
            continue
        for directory, dirs, files in os.walk(root, followlinks=False,
                onerror=lambda e: warnings.append(str(e))):
            dirs[:] = sorted(d for d in dirs if not d.startswith(".")
                             and not (Path(directory) / d).is_symlink())
            folder = Path(directory).resolve()
            if folder in seen:
                dirs[:] = []
                continue
            seen.add(folder)
            count += 1
            if count > max_dirs:
                raise ValueError("Search limit reached; use narrower --root folders")
            numbered = {}
            for filename in sorted(files):
                p = Path(filename)
                if p.suffix.lower() in (".zip", ".tar", ".gz", ".7z"):
                    archives.append(str(folder / filename))
                if p.suffix.lower() in (".tif", ".tiff") and p.stem.isdigit():
                    n = int(p.stem)
                    if 1 <= n <= 12:
                        numbered.setdefault(n, []).append(filename)
            if not numbered:
                image_files = [f for f in sorted(files) if Path(f).suffix.lower()
                               in (".tif", ".tiff", ".jpg", ".jpeg", ".nef", ".dng", ".npz")]
                if image_files:
                    other_image_folders.append({"folder": str(folder),
                        "file_count": len(image_files), "first_20_filenames": image_files[:20],
                        "status": "unrecognized_naming_no_automatic_angle_assignment"})
                continue
            missing = sorted(set(range(1, 13)) - numbered.keys())
            duplicates = {str(k): v for k, v in numbered.items() if len(v) != 1}
            complete = not missing and not duplicates
            single = set(numbered) == {1,2,3,4} and not duplicates
            entry = {"folder": str(folder), "complete_1_to_12": complete,
                     "single_1_to_4_candidate": single,
                     "layout": "three_replicates" if complete else
                               ("single_quartet_requires_confirmation" if single else "incomplete"),
                     "missing": missing, "duplicates": duplicates,
                     "files": numbered}
            candidates.append(entry)
    out = new_directory(output)
    report = {"version": VERSION, "roots": list(map(str, roots)),
              "directories_scanned": count, "warnings": warnings,
              "candidates": candidates, "other_image_folders": other_image_folders,
              "archive_files_found": archives,
              "raw_images_modified": False}
    (out / "inventory.json").write_text(json.dumps(report, indent=2) + "\n")
    for i, candidate in enumerate(candidates):
        if not (candidate["complete_1_to_12"] or candidate["single_1_to_4_candidate"]):
            continue
        groups = GROUPS if candidate["complete_1_to_12"] else [[1,2,3,4]]
        config = {
            "schema_version": 1, "site_id": None,
            "geometry_reviewed": False, "registration_reviewed": False,
            "input": {"kind": "numbered_tiffs", "directory": candidate["folder"],
                      "channel": 2,
                      "groups": [[candidate["files"][n][0] for n in g] for g in groups],
                      "rotation_deg_clockwise": [0, -45, -90, -135],
                      "shift_xy_px": [[0, 0]] * 4,
                      "stokes_method": "REVIEW_REQUIRED",
                      "comment": "Candidate angle grouping; confirm layout/order and rotation sign before use."},
            "fov": {"center_xy_px": None, "radius_60_px": None,
                    "angular_radius_deg": 60,
                    "comment": "Measure radius for 60 degrees; do not assume image-edge calibration."},
            "absorbance_parallel": 0.4, "absorbance_perpendicular": 0.04,
            "response_mode": "stokes_intensity", "invalid_policy": "error",
            "comment": "Each circular-FOV pixel contributes equally; no solid-angle weights."
        }
        (out / f"candidate_{i:04d}.json").write_text(json.dumps(config, indent=2) + "\n")
    print(json.dumps({"inventory": str(out / "inventory.json"),
                      "candidate_folders": len(candidates),
                      "complete_sets": sum(c["complete_1_to_12"] for c in candidates),
                      "single_quartet_candidates": sum(c["single_1_to_4_candidate"] for c in candidates),
                      "other_image_folders": len(other_image_folders),
                      "archive_files_found": len(archives),
                      "warnings": warnings}, indent=2))
    return report


def circle_mask(shape, center, radius):
    if len(shape) != 2 or len(center) != 2:
        raise ValueError("Expected 2D image and [x,y] centre")
    cx, cy = map(float, center)
    radius = float(radius)
    if not np.isfinite([cx, cy, radius]).all() or radius <= 0:
        raise ValueError("FOV centre/radius must be finite; radius positive")
    if cx-radius < 0 or cy-radius < 0 or cx+radius > shape[1]-1 or cy+radius > shape[0]-1:
        raise ValueError("The full 60-degree circular FOV must fit inside the image")
    y, x = np.ogrid[:shape[0], :shape[1]]
    mask = (x-cx)**2 + (y-cy)**2 <= radius**2
    if not mask.any():
        raise ValueError("Empty circular FOV")
    return mask


def align_image(image, rotation, shift):
    """Derived resampling only; positive angle is clockwise on a displayed image.

    Rotate about (width-1,height-1)/2, then translate by [dx,dy].
    Bilinear interpolation, no crop resizing, unsupported border = NaN.
    """
    from scipy.ndimage import affine_transform
    angle = math.radians(float(rotation))
    dx, dy = map(float, shift)
    if not np.isfinite([angle, dx, dy]).all():
        raise ValueError("Nonfinite registration transform")
    c, s = math.cos(angle), math.sin(angle)
    # Inverse transform in row,column coordinates for ndimage.
    matrix = np.array([[c, -s], [s, c]])
    center = (np.asarray(image.shape)-1) / 2
    offset = center - matrix @ (center + [dy, dx])
    if angle == 0 and dx == 0 and dy == 0:
        return image.copy()
    return affine_transform(image, matrix, offset, order=1,
                            mode="constant", cval=np.nan, prefilter=False)


def tiff_channel(path, channel):
    try:
        import tifffile
    except ImportError as e:
        raise ValueError("Install requirements-prc.txt for TIFF input") from e
    raw = tifffile.imread(path)
    metadata = {"path": str(path), "shape": list(raw.shape), "dtype": str(raw.dtype),
                "sha256": digest(path)}
    if raw.ndim == 2:
        if channel is not None:
            raise ValueError(f"{path}: grayscale input requires channel=null")
        image = raw
    elif raw.ndim == 3 and raw.shape[-1] in (3, 4):
        if not isinstance(channel, int) or not 0 <= channel < raw.shape[-1]:
            raise ValueError("Set a valid zero-based colour channel")
        image = raw[..., channel]
    else:
        raise ValueError(f"Unsupported TIFF layout {raw.shape}; no page/channel guessing")
    metadata["channel"] = channel
    metadata["selected_channel_max"] = float(np.nanmax(image))
    metadata["integer_max_fraction"] = (float(np.mean(image == np.iinfo(raw.dtype).max))
                                        if np.issubdtype(raw.dtype, np.integer) else None)
    return image.astype(np.float64), metadata


def inspect_inventory(inventory_path, output, preview_folder=None):
    """Read TIFF headers for the known inventory; optional first-quartet quicklook.

    No automatic analyzer order, site match or fisheye radius is assigned.
    Header inspection does not decode/test every pixel or validate TIFF integrity.
    """
    try:
        import tifffile
    except ImportError as error:
        raise ValueError("Install requirements-prc.txt for TIFF inspection") from error
    report = json.loads(Path(inventory_path).read_text())
    folders = []
    for entry in report.get("candidates", []):
        folders.append((Path(entry["folder"]),
                        [name for names in entry["files"].values() for name in names]))
    for entry in report.get("other_image_folders", []):
        # Listing the actual folder, not just the inventory's first-20 sample.
        folder = Path(entry["folder"])
        names = ([p.name for p in folder.iterdir() if p.is_file()
                  and p.suffix.lower() in (".tif", ".tiff")] if folder.is_dir() else [])
        folders.append((folder, names))
    target = Path(output).resolve()
    for folder, _ in folders:
        source = folder.resolve()
        if source == target or source in target.parents or target in source.parents:
            raise ValueError("Inspection output must be separate from all source folders")
    out = new_directory(target)
    records, errors = [], []
    for folder, names in folders:
        if not folder.is_dir():
            errors.append({"folder": str(folder), "error": "Missing image directory"})
            continue
        for name in sorted(names):
            path = folder/name
            try:
                with tifffile.TiffFile(path) as tiff:
                    page = tiff.pages[0]
                    tags = {}
                    for tag_name in ("BitsPerSample", "SamplesPerPixel", "Orientation",
                                     "ExposureTime", "FNumber", "ISOSpeedRatings", "DateTime",
                                     "Make", "Model", "Software"):
                        tag = page.tags.get(tag_name)
                        if tag is not None:
                            value = tag.value
                            # Lossless numeric/string metadata when JSON-compatible;
                            # retain a bounded representation otherwise.
                            try:
                                json.dumps(value, allow_nan=False)
                            except (TypeError, ValueError):
                                value = repr(value)[:500]
                            tags[tag_name] = value
                    records.append({"folder": str(folder), "filename": name,
                        "size_bytes": path.stat().st_size, "page_count": len(tiff.pages),
                        "page0_shape": list(page.shape), "page0_dtype": str(page.dtype),
                        "page0_photometric": str(page.photometric),
                        "page0_compression": str(page.compression), "tags": tags})
            except (OSError, ValueError, IndexError, tifffile.TiffFileError) as error:
                errors.append({"path": str(path), "error": str(error)})
    summary = {"version": VERSION, "scope": "TIFF headers only; not complete pixel integrity validation",
        "inventory_sha256": digest(inventory_path), "folder_count": len(folders),
        "tiff_headers_read": len(records), "errors": errors, "records": records,
        "raw_images_modified": False, "preview": None}
    # Preserve the header report even if a preview codec/layout is unsupported.
    (out/"tiff_inspection.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    if preview_folder is not None:
        selected = Path(preview_folder).expanduser().resolve()
        listed = {f.resolve(): names for f,names in folders}
        if selected not in listed:
            raise ValueError("Preview folder must be an exact folder in this inventory")
        # Known numbered-set first replicate, explicitly stated in output. For
        # arbitrary names show up to four only as filenames, never as angle IDs.
        names = listed[selected]
        quartet = [f"{n}.tif" for n in GROUPS[0]]
        chosen = quartet if all(n in names for n in quartet) else sorted(names)[:4]
        if not chosen:
            raise ValueError("No preview TIFFs in chosen folder")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, len(chosen), figsize=(4*len(chosen),4), squeeze=False)
        previews=[]
        for ax, name in zip(axes.flat, chosen):
            path=selected/name
            with tifffile.TiffFile(path) as tiff:
                shape=tiff.pages[0].shape
            channel=None if len(shape)==2 else 2
            image, metadata=tiff_channel(path,channel)
            stride=max(1,math.ceil(max(image.shape)/900))
            ax.imshow(image[::stride,::stride],cmap="gray")
            ax.set_title(f"{name}\n{image.shape[1]} x {image.shape[0]}, {metadata['dtype']}")
            ax.set_axis_off()
            previews.append({"filename":name,"display_channel":channel,"sha256":metadata["sha256"]})
        fig.suptitle("UNALIGNED INPUT QUICKLOOK — individual display scales; no FOV calibration")
        fig.tight_layout()
        fig.savefig(out/"input_quicklook.png",dpi=130)
        plt.close(fig)
        summary["preview"]={"folder":str(selected),"images":previews,
            "label":"unaligned inputs; not analyzer-order or registration validation"}
    (out/"tiff_inspection.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"output":str(out),"tiff_headers_read":len(records),
                      "errors":errors,"preview_created":summary["preview"] is not None},indent=2))
    return summary


def combine_replicates(replicates, method):
    """Replicates are already co-registered [I0,I45,I90,I135] arrays.

    legacy_dop_aop reproduces the R DoP display masks and arithmetic AoP
    averaging, then reconstructs Q/U; it is not evidence that this intermediate
    was used in the historical spider integration. linear_stokes is explicit
    alternate processing. Both are kept separate in metadata.
    """
    if method not in ("legacy_dop_aop", "linear_stokes"):
        raise ValueError("Choose stokes_method: legacy_dop_aop or linear_stokes")
    si = sq = su = sd = sa = None
    qa = []
    shape = None
    for quartet in replicates:
        if len(quartet) != 4 or any(x.ndim != 2 for x in quartet):
            raise ValueError("Four 2D images required per replicate")
        if shape is None:
            shape = quartet[0].shape
        if any(x.shape != shape for x in quartet):
            raise ValueError("Analyzer/replicate image dimensions disagree")
        z, f, n, o = quartet
        i, q, u = z+n, z-n, f-o
        with np.errstate(divide="ignore", invalid="ignore"):
            qn, un = (z-n)/(z+n), (f-o)/(f+o)
            dop, aop = np.hypot(qn, un), 0.5*np.arctan2(un, qn)
            closure = ((z+n)-(f+o)) / ((z+n)+(f+o))
        qa.append({"nonfinite_dop_pixels": int(np.sum(~np.isfinite(dop))),
                   "dop_gt_0p95_pixels": int(np.sum(dop > .95)),
                   "closure": closure})
        if method == "legacy_dop_aop":
            # Preserve zeroing from supplied R; keep NaN angles/unsupported I.
            dop = np.where(np.isnan(dop) | (dop > .95), 0, dop)
        if si is None:
            si, sq, su, sd, sa = [np.zeros(shape) for _ in range(5)]
        si += i
        sq += q
        su += u
        sd += dop
        sa += aop
    count = len(qa)
    if not count:
        raise ValueError("No replicate images")
    si /= count
    if method == "legacy_dop_aop":
        d, a = sd/count, sa/count
        sq, su = si*d*np.cos(2*a), si*d*np.sin(2*a)
    else:
        sq /= count
        su /= count
    return si, sq, su, qa


def load_sky(config, base):
    inp = config["input"]
    kind = inp["kind"]
    records, qa = [], []
    if kind in ("stokes_npz", "dop_aop_npz"):
        path = relative_path(inp["path"], base)
        records.append({"path": str(path), "sha256": digest(path)})
        with np.load(path, allow_pickle=False) as data:
            if kind == "stokes_npz":
                i, q, u = [np.asarray(data[k], dtype=np.float64) for k in ("I", "Q", "U")]
            else:
                d = np.asarray(data["DoP"], dtype=np.float64)
                a = np.asarray(data["AoP"], dtype=np.float64)
                if d.ndim != 2 or a.shape != d.shape:
                    raise ValueError("DoP/AoP must be identically shaped 2D arrays")
                units = inp.get("aop_unit")
                if units == "degrees":
                    a = np.deg2rad(a)
                elif units != "radians":
                    raise ValueError("Specify aop_unit radians or degrees")
                if "I" in data:
                    i = np.asarray(data["I"], dtype=np.float64)
                    if i.shape != d.shape:
                        raise ValueError("I and DoP/AoP dimensions disagree")
                elif config["response_mode"] == "normalized_direction":
                    i = np.ones_like(d)
                else:
                    raise ValueError("Intensity-mode DoP/AoP input must include I")
                q, u = i*d*np.cos(2*a), i*d*np.sin(2*a)
    elif kind in ("numbered_tiffs", "tiff_groups"):
        directory = relative_path(inp["directory"], base)
        groups = inp["groups"]
        if not 1 <= len(groups) <= 3 or any(len(g) != 4 for g in groups):
            raise ValueError("TIFF mode requires one to three explicit four-angle groups")
        paths = [relative_path(f, directory) for g in groups for f in g]
        if len(set(paths)) != 4*len(groups):
            raise ValueError("All TIFF inputs must be distinct; do not duplicate a replicate")
        rotations, shifts = inp["rotation_deg_clockwise"], inp["shift_xy_px"]
        if len(rotations) != 4 or len(shifts) != 4:
            raise ValueError("One rotation/translation per analyzer required")
        def generator():
            for group in groups:
                quartet = []
                for index, filename in enumerate(group):
                    path = relative_path(filename, directory)
                    image, record = tiff_channel(path, inp["channel"])
                    records.append(record)
                    quartet.append(align_image(image, rotations[index], shifts[index]))
                yield quartet
        i, q, u, qa = combine_replicates(generator(), inp["stokes_method"])
    else:
        raise ValueError(f"Unsupported input kind: {kind}")
    if i.ndim != 2 or i.shape != q.shape or i.shape != u.shape:
        raise ValueError("I/Q/U must be identically shaped 2D arrays")
    return i, q, u, records, qa


def eye_result(i, q, u, mask, a_parallel=.4, a_perpendicular=.04,
               response_mode="stokes_intensity", invalid_policy="error"):
    if not np.isfinite([a_parallel, a_perpendicular]).all() or not 0 <= a_perpendicular <= a_parallel:
        raise ValueError("Require finite 0 <= perpendicular <= parallel absorbance")
    physical = (np.isfinite(i) & np.isfinite(q) & np.isfinite(u) & (i > 0)
                & (np.hypot(q, u) <= i*(1+1e-10)))
    invalid = mask & ~physical
    if invalid_policy not in ("error", "exclude"):
        raise ValueError("invalid_policy must be error or exclude")
    if invalid.any() and invalid_policy == "error":
        raise ValueError(f"{int(invalid.sum())} invalid FOV pixels; inspect before choosing exclusion")
    selected = mask & physical
    count = int(selected.sum())
    if not count:
        raise ValueError("No valid FOV pixels")
    iv, qv, uv = [x[selected] for x in (i, q, u)]
    if response_mode == "normalized_direction":
        qv, uv = qv/iv, uv/iv
        iv = np.ones_like(iv)
    elif response_mode != "stokes_intensity":
        raise ValueError("Choose stokes_intensity or normalized_direction response_mode")
    total_i, total_q, total_u = map(float, (iv.sum(), qv.sum(), uv.sum()))
    r, t = 10**(-a_parallel), 10**(-a_perpendicular)
    b, a = 1-(r+t)/2, (t-r)/2
    mean = b*total_i
    amplitude = a*math.hypot(total_q, total_u)
    maximum, minimum = mean+amplitude, mean-amplitude
    if minimum <= 0 or not np.isfinite([maximum, minimum]).all():
        raise ValueError("Nonpositive/nonfinite integrated eye response")
    angle = (0.5*math.atan2(-total_u, -total_q)) % math.pi if amplitude > 0 else None
    return {"prc": math.log(maximum/minimum), "response_max_sum": maximum,
            "response_min_sum": minimum, "response_max_mean": maximum/count,
            "response_min_mean": minimum/count, "valid_pixels": count,
            "fov_pixels": int(mask.sum()), "excluded_pixels": int(invalid.sum()),
            "excluded_fraction": float(invalid.sum()/mask.sum()),
            "sum_I": total_i, "sum_Q": total_q, "sum_U": total_u,
            "effective_integrated_DoLP": math.hypot(total_q,total_u)/total_i,
            "max_axis_deg_mod_180": None if angle is None else math.degrees(angle),
            "min_axis_deg_mod_180": None if angle is None else (math.degrees(angle)+90)%180,
            "axis_convention": "Q cos(2 theta)+U sin(2 theta); zero along +Q reference",
            "absorbance_parallel": a_parallel, "absorbance_perpendicular": a_perpendicular,
            "response_mode": response_mode, "spatial_weights": "equal included pixels"}


def reference_comparison(site_id, prc):
    if site_id is None:
        return {"status": "no_site_id_selected"}
    with (ROOT / "data/reference/prc_spider_18_sites_verified.csv").open() as stream:
        rows = [row for row in csv.DictReader(stream) if row["site_id"] == site_id]
    if len(rows) != 1:
        raise ValueError("site_id must uniquely match verified spider reference")
    expected = float(rows[0]["prc_sigg"])
    return {"site_id": site_id, "location": rows[0]["location"],
            "reference_rounded_prc": expected, "difference": prc-expected,
            "within_3_decimal_rounding": abs(prc-expected) <= .000500000001,
            "status": "numeric_comparison_only_not_validation_of_geometry_or_image_identity"}


def preview_plot(out, i, q, u, mask, result, qa):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    stride = max(1, math.ceil(max(i.shape)/900))
    with np.errstate(divide="ignore", invalid="ignore"):
        dop, aop = np.hypot(q,u)/i, np.rad2deg(.5*np.arctan2(u,q))
    arrays = [i, dop, aop]
    titles = ["Derived intensity", "DoLP (unclipped values retained)", "AoP, degrees"]
    if qa:
        arrays += [x["closure"] for x in qa]
        titles += [f"Replicate {j+1}: analyzer-pair closure" for j in range(len(qa))]
    rows = 2 if qa else 1
    fig, axes = plt.subplots(rows, 3, figsize=(13, 4*rows), squeeze=False)
    for j, ax in enumerate(axes.flat):
        if j < len(arrays):
            artist = ax.imshow(arrays[j][::stride,::stride], cmap="twilight" if j == 2 else "viridis")
            ax.contour(mask[::stride,::stride], levels=[.5], colors="red", linewidths=.8)
            ax.set_title(titles[j])
            fig.colorbar(artist, ax=ax, shrink=.7)
        ax.set_axis_off()
    fig.suptitle(f"RECONSTRUCTION QC — {result['status']} — PRC {result['prc']:.6f}")
    fig.tight_layout()
    fig.savefig(out / "qc_preview.png", dpi=140)
    plt.close(fig)


def run(manifest, output, preview_only=False):
    manifest = Path(manifest).resolve()
    config = json.loads(manifest.read_text())
    if config.get("schema_version") != 1:
        raise ValueError("Unsupported manifest schema")
    if not preview_only and not (config.get("geometry_reviewed") is True
                                and config.get("registration_reviewed") is True):
        raise ValueError("Review geometry and registration first; use --preview-only for QC")
    fov = config["fov"]
    if fov.get("angular_radius_deg") != 60:
        raise ValueError("This reconstruction requires the specified 60-degree angular radius")
    if fov.get("center_xy_px") is None or fov.get("radius_60_px") is None:
        raise ValueError("Set measured center_xy_px and radius_60_px; no inferred fisheye calibration")
    # Reject source/output overlap before creating outputs.
    inp = config["input"]
    source = relative_path(inp.get("directory", inp.get("path", ".")), manifest.parent)
    target = Path(output).resolve()
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {target}")
    if source == target or source in target.parents or target in source.parents:
        raise ValueError("Output must be separate from input image directory/file")
    i, q, u, records, qa = load_sky(config, manifest.parent)
    mask = circle_mask(i.shape, fov["center_xy_px"], fov["radius_60_px"])
    result = eye_result(i, q, u, mask, config["absorbance_parallel"],
                        config["absorbance_perpendicular"], config["response_mode"],
                        config["invalid_policy"])
    result["status"] = "UNREVIEWED_PREVIEW" if preview_only else "RECONSTRUCTED_NOT_HISTORICALLY_VALIDATED"
    result["replicate_count"] = len(qa) if inp["kind"] in ("numbered_tiffs", "tiff_groups") else None
    result["reference"] = reference_comparison(config.get("site_id"), result["prc"])
    for record in qa:
        values = record["closure"][mask]
        values = values[np.isfinite(values)]
        record["closure_fov_abs_median"] = float(np.median(abs(values))) if len(values) else None
        record["closure_fov_abs_p95"] = float(np.quantile(abs(values), .95)) if len(values) else None
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
    dependencies = {}
    for package in ("numpy", "scipy", "matplotlib", "tifffile"):
        try:
            dependencies[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            dependencies[package] = None
    provenance = {"version": VERSION, "created_utc": datetime.now(timezone.utc).isoformat(),
                  "git_commit": git.stdout.strip() or None,
                  "script_sha256": digest(__file__), "manifest_sha256": digest(manifest),
                  "config": config, "inputs": records,
                  "qa": [{k:v for k,v in rec.items() if k != "closure"} for rec in qa],
                  "python_version": platform.python_version(), "dependencies": dependencies,
                  "raw_images_modified": False,
                  "limitations": ["60-degree pixel radius must be independently calibrated",
                    "global eye-axis optimisation is not exposure image registration",
                    "common circular FOV and global image-plane axes; no inferred eye curvature",
                    "no spectral correction, no VIIRS refit",
                    "historical validation requires verified matching original images and output"]}
    # Refuse overwrite; raw inputs are opened read-only, derived arrays float64.
    out = new_directory(target)
    np.savez_compressed(out / "derived_stokes.npz", I=i, Q=q, U=u, fov_mask=mask)
    (out / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False)+"\n")
    b = (result["response_max_sum"]+result["response_min_sum"])/2
    amp = (result["response_max_sum"]-result["response_min_sum"])/2
    angle = math.radians(result["max_axis_deg_mod_180"] or 0)
    with (out / "eye_orientation_curve.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["axis_deg", "eye_1_sum", "eye_2_sum", "signed_log_ratio"])
        for theta_deg in np.arange(0, 180, .5):
            modulation = amp*math.cos(2*(math.radians(theta_deg)-angle))
            first, second = b+modulation, b-modulation
            writer.writerow([theta_deg, first, second, math.log(first/second)])
    preview_plot(out, i, q, u, mask, result, qa)
    print(json.dumps({"output": str(out), **result}, indent=2, allow_nan=False))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    inventory = commands.add_parser("discover", help="Locate numbered 1–12 TIFF sets")
    inventory.add_argument("--root", type=Path, action="append", required=True)
    inventory.add_argument("--output", type=Path, required=True)
    inventory.add_argument("--max-dirs", type=int, default=100000)
    inspect = commands.add_parser("inspect", help="Audit TIFF headers and optionally preview one folder")
    inspect.add_argument("--inventory", type=Path, required=True)
    inspect.add_argument("--output", type=Path, required=True)
    inspect.add_argument("--preview-folder", type=Path)
    process = commands.add_parser("run", help="Process one explicitly selected manifest")
    process.add_argument("--manifest", type=Path, required=True)
    process.add_argument("--output", type=Path, required=True)
    process.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "discover":
            discover(args.root, args.output, args.max_dirs)
        elif args.command == "inspect":
            inspect_inventory(args.inventory, args.output, args.preview_folder)
        else:
            run(args.manifest, args.output, args.preview_only)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, f"Error: {error}\n")


if __name__ == "__main__":
    main()
