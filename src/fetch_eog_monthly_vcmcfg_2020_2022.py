#!/usr/bin/env python3
"""Download and extract monthly EOG VCMCFG for 2020-2022.

This stage fills the missing VCMCFG half of the fixed-NUTS2016 monthly
VCMCFG/VCMSLCFG product-atlas diagnostic.

The existing July 2020/2021/2022 VCMCFG pilot files are reused without
re-downloading them. All other missing months are fetched from the official
EOG monthly VNL V1 repository.

Authentication follows the same established EOG OAuth/OpenID workflow already
used by fetch_eog_monthly_vcmslcfg_2021.py. Credentials are used only in memory
and are never written to disk or committed.
"""

from __future__ import annotations

import argparse
import csv
from html.parser import HTMLParser
from pathlib import Path
import re
import shutil
import tarfile
from urllib.parse import urljoin

import rasterio

from fetch_eog_monthly_vcmslcfg_2021 import (
    EOGAuth,
    download_archive,
    read_url_text,
)
from project_config import PROJECT_ROOT, load_config


PRODUCT = "vcmcfg"
ARCHIVE_RE = re.compile(
    r"^SVDNB_(?P<satellite>npp|j01)_(?P<start>\d{8})-(?P<end>\d{8})_"
    r"(?P<tile>[^_]+)_vcmcfg_v10_.*\.tgz$",
    re.IGNORECASE,
)


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def find_archive_in_listing(
    listing_url: str,
    html: str,
    *,
    year: int,
    month: int,
    tile: str,
) -> tuple[str, str]:
    parser = LinkParser()
    parser.feed(html)
    prefix = f"SVDNB_"
    ym = f"_{year}{month:02d}"
    candidates: list[tuple[str, str]] = []
    for href in parser.hrefs:
        name = Path(href.split("?", 1)[0]).name
        lower = name.lower()
        if (
            name.startswith(prefix)
            and ym in name
            and f"_{tile}_{PRODUCT}_v10_" in lower
            and lower.endswith(".tgz")
        ):
            candidates.append((name, urljoin(listing_url, href)))

    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected exactly one {tile} {PRODUCT.upper()} archive in "
            f"{listing_url}; found {len(candidates)}: {[x[0] for x in candidates]}"
        )
    return candidates[0]


def determine_access_mode(
    listing_url: str,
    auth: EOGAuth,
    *,
    year: int,
    month: int,
    tile: str,
):
    try:
        html = read_url_text(listing_url, auth, bearer=False)
        archive = find_archive_in_listing(
            listing_url, html, year=year, month=month, tile=tile
        )
        print("EOG repository permits direct programmatic access.", flush=True)
        return False, archive
    except Exception as direct_exc:
        print(
            "Direct repository access did not expose the requested archive; "
            "trying official OAuth/OpenID access.",
            flush=True,
        )
        try:
            auth.obtain_token()
            html = read_url_text(listing_url, auth, bearer=True)
            archive = find_archive_in_listing(
                listing_url, html, year=year, month=month, tile=tile
            )
            print("EOG OAuth/OpenID access: OK", flush=True)
            return True, archive
        except Exception as oauth_exc:
            raise RuntimeError(
                "Unable to access the EOG repository programmatically.\n"
                f"Direct access error: {direct_exc}\n"
                f"OAuth access error: {oauth_exc}"
            ) from oauth_exc


def canonical_root(template: str, year: int) -> Path:
    return PROJECT_ROOT / template.format(year=year)


def pilot_files(year: int, tile: str, layers: list[str]) -> dict[str, Path]:
    pilot = PROJECT_ROOT / "data/raw/viirs/monthly_july_pilot"
    if not pilot.exists():
        return {}
    result: dict[str, Path] = {}
    for layer in layers:
        matches = sorted(
            pilot.glob(
                f"SVDNB_*_{year}07*-*_{tile}_{PRODUCT}_v10_*.{layer}.tif"
            )
        )
        if len(matches) == 1:
            result[layer] = matches[0]
    return result if len(result) == len(layers) else {}


def link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.stat().st_size > 0:
        return
    try:
        dst.hardlink_to(src)
        print(f"  hard-linked existing July pilot file: {dst.name}", flush=True)
    except OSError:
        shutil.copy2(src, dst)
        print(f"  copied existing July pilot file: {dst.name}", flush=True)


def local_archive_candidates(
    year: int,
    month: int,
    tile: str,
    archive_dir: Path,
) -> list[Path]:
    patterns = [
        f"SVDNB_npp_{year}{month:02d}*_{tile}_{PRODUCT}_v10_*.tgz",
        f"SVDNB_j01_{year}{month:02d}*_{tile}_{PRODUCT}_v10_*.tgz",
    ]
    candidates: list[Path] = []
    for root in [archive_dir, Path.home() / "Downloads"]:
        if not root.exists():
            continue
        for pattern in patterns:
            candidates.extend(sorted(root.glob(pattern)))
    unique: list[Path] = []
    seen = set()
    for p in candidates:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique.append(p)
    return unique


def extract_layers(
    archive: Path,
    tif_dir: Path,
    layers: list[str],
    tile: str,
) -> dict[str, Path]:
    tif_dir.mkdir(parents=True, exist_ok=True)
    extracted: dict[str, Path] = {}
    with tarfile.open(archive, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        for layer in layers:
            matches = [
                m for m in members
                if Path(m.name).name.endswith(f".{layer}.tif")
                and f"_{tile}_{PRODUCT}_v10_" in Path(m.name).name.lower()
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"Expected one {layer} layer in {archive.name}; "
                    f"found {len(matches)}."
                )
            member = matches[0]
            name = Path(member.name).name
            dst = tif_dir / name
            if dst.exists() and dst.stat().st_size > 0:
                extracted[layer] = dst
                continue
            part = dst.with_suffix(dst.suffix + ".part")
            if part.exists():
                part.unlink()
            source = tf.extractfile(member)
            if source is None:
                raise RuntimeError(f"Cannot read {member.name} from {archive}")
            with source, part.open("wb") as out:
                shutil.copyfileobj(source, out, length=8 * 1024 * 1024)
            part.replace(dst)
            extracted[layer] = dst
            print(f"  extracted {name}", flush=True)
    return extracted


def validate_extracted(
    paths: dict[str, Path],
    expected_crs: str,
) -> dict[str, dict]:
    ref = None
    metadata: dict[str, dict] = {}
    for layer, path in sorted(paths.items()):
        with rasterio.open(path) as ds:
            if ds.count != 1:
                raise RuntimeError(f"Expected one band: {path}")
            if ds.crs is None or ds.crs.to_string() != expected_crs:
                raise RuntimeError(
                    f"Unexpected CRS for {path.name}: {ds.crs}; "
                    f"expected {expected_crs}"
                )
            sig = (
                ds.width,
                ds.height,
                tuple(ds.transform),
                tuple(ds.bounds),
                ds.crs.to_string(),
            )
            if ref is None:
                ref = sig
            elif sig != ref:
                raise RuntimeError(
                    f"Extracted layers are not aligned for {path.name}."
                )
            metadata[layer] = {
                "filename": path.name,
                "size_bytes": int(path.stat().st_size),
                "width": ds.width,
                "height": ds.height,
                "crs": ds.crs.to_string(),
                "nodata": ds.nodata,
            }
    return metadata


def parse_archive_dates(name: str) -> tuple[str, str, str]:
    m = ARCHIVE_RE.match(name)
    if not m:
        raise RuntimeError(f"Unexpected archive filename: {name}")
    return m.group("start"), m.group("end"), m.group("satellite").lower()


def write_manifest(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "year", "month", "start_date", "end_date", "satellite",
        "archive_name", "archive_source",
        "avg_rade9h_file", "cf_cvg_file", "cvg_file",
        "width", "height", "crs", "validated",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def existing_layers(
    tif_dir: Path,
    *,
    year: int,
    month: int,
    tile: str,
    layers: list[str],
) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for layer in layers:
        matches = sorted(
            tif_dir.glob(
                f"SVDNB_*_{year}{month:02d}*-*_{tile}_{PRODUCT}_v10_*.{layer}.tif"
            )
        )
        if len(matches) == 1:
            result[layer] = matches[0]
        elif len(matches) > 1:
            raise RuntimeError(
                f"Multiple canonical {PRODUCT} {year}-{month:02d} {layer} files: "
                f"{[x.name for x in matches]}"
            )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["downloads"]["eog_monthly_vcmcfg_2020_2022"]

    years = [int(x) for x in cfg["years"]]
    months = [int(x) for x in cfg["months"]]
    tile = str(cfg["tile"])
    layers = [str(x) for x in cfg["layers"]]
    base_url = str(cfg["base_url"]).rstrip("/")
    raw_dir_template = str(cfg["output_dir_template"])
    keep_archives = bool(cfg.get("keep_archives", False))

    required_layers = {"avg_rade9h", "cf_cvg", "cvg"}
    if set(layers) != required_layers:
        raise ValueError(
            "VCMCFG downloader expects exactly "
            f"{sorted(required_layers)}; configured {sorted(layers)}."
        )

    auth = EOGAuth()
    access_mode_known = False
    bearer = False
    rows: list[dict] = []

    for year in years:
        out_root = canonical_root(raw_dir_template, year)
        archive_dir = out_root / "archives"
        tif_dir = out_root / "tif"
        manifest_path = out_root / "manifest.csv"
        archive_dir.mkdir(parents=True, exist_ok=True)
        tif_dir.mkdir(parents=True, exist_ok=True)

        # Reuse the validated July pilot files instead of downloading July again.
        july = pilot_files(year, tile, layers)
        if july:
            for layer, src in july.items():
                link_or_copy(src, tif_dir / src.name)

        year_rows: list[dict] = []
        for month in months:
            print(
                f"\n===== {year}-{month:02d} {PRODUCT.upper()} {tile} =====",
                flush=True,
            )

            existing = existing_layers(
                tif_dir,
                year=year,
                month=month,
                tile=tile,
                layers=layers,
            )
            if len(existing) == len(layers):
                meta = validate_extracted(existing, config["nuts"]["crs"])
                archive_name = re.sub(
                    r"\.avg_rade9h\.tif$",
                    ".tgz",
                    existing["avg_rade9h"].name,
                )
                start, end, sat = parse_archive_dates(archive_name)
                row = {
                    "year": year,
                    "month": month,
                    "start_date": start,
                    "end_date": end,
                    "satellite": sat,
                    "archive_name": archive_name,
                    "archive_source": (
                        "existing_july_pilot" if month == 7 else "existing_extracted"
                    ),
                    "avg_rade9h_file": existing["avg_rade9h"].name,
                    "cf_cvg_file": existing["cf_cvg"].name,
                    "cvg_file": existing["cvg"].name,
                    "width": meta["avg_rade9h"]["width"],
                    "height": meta["avg_rade9h"]["height"],
                    "crs": meta["avg_rade9h"]["crs"],
                    "validated": True,
                }
                year_rows.append(row)
                write_manifest(manifest_path, year_rows)
                print("  already extracted and validated; skipping download", flush=True)
                continue

            # Prefer archives already downloaded manually.  This keeps the
            # established browser-download workflow independent of EOG
            # programmatic credentials and avoids any network/authentication
            # attempt when the required .tgz is already in ~/Downloads.
            local = local_archive_candidates(year, month, tile, archive_dir)
            if len(local) > 1:
                raise RuntimeError(
                    f"More than one local archive candidate for {year}-{month:02d}: "
                    f"{[str(x) for x in local]}"
                )

            if local:
                archive = local[0]
                archive_source = "local_existing"
                print(f"  reusing local archive: {archive}", flush=True)
            else:
                listing_url = f"{base_url}/{year}/{year}{month:02d}/{PRODUCT}/"
                if not access_mode_known:
                    bearer, first_archive = determine_access_mode(
                        listing_url,
                        auth,
                        year=year,
                        month=month,
                        tile=tile,
                    )
                    access_mode_known = True
                else:
                    html = read_url_text(listing_url, auth, bearer=bearer)
                    first_archive = find_archive_in_listing(
                        listing_url,
                        html,
                        year=year,
                        month=month,
                        tile=tile,
                    )

                archive_name, archive_url = first_archive
                archive = archive_dir / archive_name
                if not archive.exists():
                    download_archive(
                        archive_url,
                        archive,
                        auth,
                        bearer=bearer,
                    )
                archive_source = "downloaded"

            extracted = extract_layers(archive, tif_dir, layers, tile)
            meta = validate_extracted(extracted, config["nuts"]["crs"])
            start, end, sat = parse_archive_dates(archive.name)

            row = {
                "year": year,
                "month": month,
                "start_date": start,
                "end_date": end,
                "satellite": sat,
                "archive_name": archive.name,
                "archive_source": archive_source,
                "avg_rade9h_file": extracted["avg_rade9h"].name,
                "cf_cvg_file": extracted["cf_cvg"].name,
                "cvg_file": extracted["cvg"].name,
                "width": meta["avg_rade9h"]["width"],
                "height": meta["avg_rade9h"]["height"],
                "crs": meta["avg_rade9h"]["crs"],
                "validated": True,
            }
            year_rows.append(row)
            write_manifest(manifest_path, year_rows)

            if (
                archive_source == "downloaded"
                and not keep_archives
                and archive.parent.resolve() == archive_dir.resolve()
            ):
                archive.unlink()
                print(
                    "  removed downloaded archive after validated extraction",
                    flush=True,
                )

        if len(year_rows) != len(months):
            raise RuntimeError(
                f"{year}: expected {len(months)} completed months; "
                f"manifest has {len(year_rows)}."
            )
        rows.extend(year_rows)

    print("\nEOG MONTHLY VCMCFG DOWNLOAD COMPLETE")
    print(f"Years: {years}")
    print(f"Completed product-months: {len(rows)}")
    print(f"Reused July pilot months: {len(years)}")
    print(f"Downloaded/validated missing months: {len(rows) - len(years)}")
    print("Run the monthly atlas inventory next and require 36/36 for both products.")


if __name__ == "__main__":
    main()
