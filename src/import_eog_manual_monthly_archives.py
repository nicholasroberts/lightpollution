#!/usr/bin/env python3
"""Import already-downloaded monthly EOG VNL V1 .tgz archives from ~/Downloads.

Read only from archives, extract just avg_rade9h and cf_cvg, validate matched
GeoTIFF grids, and write separate manual import manifests. Source archives
are never deleted or moved. This is NOT an automated network downloader.

Examples:
    python src/import_eog_manual_monthly_archives.py --dry-run
    python src/import_eog_manual_monthly_archives.py
    python src/import_eog_manual_monthly_archives.py --years 2019 2023 --product vcmslcfg
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import tarfile
from collections import defaultdict
from pathlib import Path

import rasterio

from project_config import PROJECT_ROOT

ARCHIVE = re.compile(
    r"^SVDNB_(?:npp|j01)_(?P<start>\d{8})-(?P<end>\d{8})_"
    r"(?P<tile>\d{2}N\d{3}[WE])_(?P<product>vcmslcfg|vcmcfg)_v10_.*\.tgz$",
    re.IGNORECASE,
)
LAYERS = ("avg_rade9h", "cf_cvg")
FIELDNAMES = [
    "product", "year", "month", "start_date", "end_date",
    "archive_name", "archive_source", "avg_rade9h_file",
    "cf_cvg_file", "validated",
]


def candidates(root, tile, product, years):
    out = defaultdict(list)
    for archive in sorted(root.glob("*.tgz")):
        m = ARCHIVE.match(archive.name)
        if m is None:
            continue
        year = int(m["start"][:4])
        mo = int(m["start"][4:6])
        if (m["tile"].upper() == tile and m["product"].lower() == product
                and year in years and 1 <= mo <= 12):
            out[(year, mo)].append((archive, m))
    return out


def extract(archive, outdir, tile, product, dry_run):
    """Never extract arbitrary TAR paths: copy named expected members to basename."""
    with tarfile.open(archive, "r:gz") as tf:
        found = {}
        for member in tf.getmembers():
            if not member.isfile():
                continue
            name = Path(member.name).name
            if f"_{tile}_{product}_v10_" not in name.lower():
                continue
            for layer in LAYERS:
                if name.lower().endswith(f".{layer}.tif"):
                    if layer in found:
                        raise RuntimeError(f"Duplicate {layer} in {archive.name}")
                    found[layer] = member
        if set(found) != set(LAYERS):
            raise RuntimeError(f"Missing required layers in {archive.name}: "
                               f"found {sorted(found)}")
        result = {}
        for layer in LAYERS:
            member = found[layer]
            dst = outdir / Path(member.name).name
            result[layer] = dst
            if dry_run:
                print(f"    would extract: {dst.relative_to(PROJECT_ROOT)}")
                continue
            if dst.exists() and dst.stat().st_size > 0:
                print(f"    existing: {dst.name}")
                continue
            outdir.mkdir(parents=True, exist_ok=True)
            tmp = dst.with_name(dst.name + ".part")
            src = tf.extractfile(member)
            if src is None:
                raise RuntimeError(f"Cannot read TAR member {member.name}")
            try:
                with src, tmp.open("wb") as w:
                    shutil.copyfileobj(src, w, length=8 * 1024 * 1024)
                tmp.replace(dst)
            except Exception:
                tmp.unlink(missing_ok=True)
                raise
            print(f"    extracted: {dst.name}")
    return result


def validate(layers):
    signatures = []
    for name in LAYERS:
        with rasterio.open(layers[name]) as ds:
            if ds.count != 1 or ds.crs is None:
                raise ValueError(f"Invalid 1-band georeferenced raster: {layers[name]}")
            signatures.append(
                (ds.width, ds.height, ds.crs.to_string(), tuple(ds.transform)))
    if signatures[0] != signatures[1]:
        raise ValueError("Radiance and cloud-free rasters do not share a grid")


def append_manifest(path, row):
    old = []
    if path.exists():
        with path.open(newline="", encoding="utf-8") as f:
            old = list(csv.DictReader(f))
    # Repeat executions are idempotent for each product-year-month.
    old = [r for r in old if not (
        r.get("product") == row["product"] and
        r.get("year") == row["year"] and
        r.get("month") == row["month"])]
    old.append(row)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".csv.part")
    with temp.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(old)
    temp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive-dir", type=Path,
                   default=Path.home() / "Downloads")
    p.add_argument("--tile", default="75N060W")
    p.add_argument("--product", choices=("vcmslcfg", "vcmcfg"),
                   default="vcmslcfg")
    p.add_argument("--years", nargs="+", type=int,
                   default=[2015, 2016, 2017, 2018, 2019, 2023, 2024])
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    if not args.archive_dir.is_dir():
        p.error(f"Archive directory missing: {args.archive_dir}")
    tile = args.tile.upper()
    product = args.product.lower()
    items = candidates(args.archive_dir, tile, product, set(args.years))
    if not items:
        print("No matching .tgz archives found. Nothing to import.")
        return
    processed = 0
    for (year, mo), matches in sorted(items.items()):
        if len(matches) != 1:
            raise RuntimeError(
                f"Ambiguous {year}-{mo:02d} archives: "
                f"{[x[0].name for x in matches]}")
        archive, stamp = matches[0]
        outdir = (PROJECT_ROOT / "data/raw/viirs" /
                  f"monthly_{product}_{year}" / "tif")
        print(f"\n{year}-{mo:02d} {product}: {archive.name}")
        layers = extract(archive, outdir, tile, product, args.dry_run)
        if not args.dry_run:
            validate(layers)
            row = {
                "product": product, "year": str(year), "month": str(mo),
                "start_date": stamp["start"], "end_date": stamp["end"],
                "archive_name": archive.name, "archive_source": "manual_import",
                "avg_rade9h_file": layers["avg_rade9h"].name,
                "cf_cvg_file": layers["cf_cvg"].name, "validated": "True"}
            append_manifest(outdir.parent / "manual_import_manifest.csv", row)
            print("    grid validated")
        processed += 1
    print(f"\n{'Would process' if args.dry_run else 'Processed'} "
          f"{processed} {product} month archive(s).")
    print("Archives in Downloads retained; annual products unchanged.")


if __name__ == "__main__":
    main()
