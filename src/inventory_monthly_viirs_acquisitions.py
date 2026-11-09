#!/usr/bin/env python3
"""Read-only EcoVis EOG monthly inventory, distinguishing observation dates from build dates.

Run: python src/inventory_monthly_viirs_acquisitions.py
No downloading, unzipping or modification to source products. Prints a compact
year/product matrix and summarizes downloader manifest provenance when present.
Only the generated inventory CSV is written to data/processed/.
"""
from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from project_config import PROJECT_ROOT

ROOT = PROJECT_ROOT
RAW = ROOT / "data/raw/viirs"
OUT = ROOT / "data/processed/local_lighting_interventions/monthly_acquisition_inventory.csv"
TILE = "75N060W"
YEARS = range(2015, 2025)
PRODUCTS = ("vcmslcfg", "vcmcfg")
# Crucial: first YYYYMMDD after satellite, NEVER an EOG _cYYYYMMDD creation stamp.
NAME = re.compile(
    r"^SVDNB_(?:npp|j01)_(?P<start>20\d{6})-(?P<end>20\d{6})_"
    r"(?P<tile>\d{2}N\d{3}[WE])_(?P<product>vcmslcfg|vcmcfg)_v\d+_",
    re.IGNORECASE,
)


def source_identity(path):
    m = NAME.match(path.name)
    if not m or m.group("tile").upper() != TILE:
        return None
    yr, mo = int(m.group("start")[:4]), int(m.group("start")[4:6])
    if yr not in YEARS or not 1 <= mo <= 12:
        return None
    return m.group("product").lower(), yr, mo


def main():
    layers = defaultdict(lambda: {"rad": [], "cf": []})
    archives = defaultdict(list)
    manifests = []
    ignored_names = 0
    if not RAW.exists():
        raise FileNotFoundError(RAW)

    for p in RAW.rglob("*"):
        if not p.is_file():
            continue
        if p.name.lower().endswith("manifest.csv"):
            manifests.append(p)
        suffix = p.suffix.lower()
        if suffix not in {".tif", ".tiff", ".tgz"}:
            continue
        key = source_identity(p)
        if key is None:
            ignored_names += 1
            continue
        if suffix == ".tgz":
            archives[key].append(str(p))
        elif "avg_rade9" in p.name.lower():
            layers[key]["rad"].append(str(p))
        elif "cf_cvg" in p.name.lower():
            layers[key]["cf"].append(str(p))

    downloads = Path.home() / "Downloads"
    if downloads.exists():
        for p in downloads.glob("*.tgz"):
            key = source_identity(p)
            if key:
                archives[key].append(str(p))

    rows = []
    for year in YEARS:
        for product in PRODUCTS:
            complete = 0
            missing = []
            pending = []
            for mo in range(1, 13):
                key = (product, year, mo)
                d = layers[key]
                rad, cf = bool(d["rad"]), bool(d["cf"])
                ready = rad and cf
                complete += int(ready)
                if not ready:
                    missing.append(f"{mo:02d}")
                    if archives[key]:
                        pending.append(f"{mo:02d}")
                rows.append({
                    "year": year, "month": mo, "product": product,
                    "ready": int(ready),
                    "radiance_files": len(d["rad"]),
                    "cloudfree_files": len(d["cf"]),
                    "archives_found": len(archives[key]),
                    "radiance_path": d["rad"][0] if rad else "",
                    "cloudfree_path": d["cf"][0] if cf else "",
                    "archive_path": archives[key][0] if archives[key] else "",
                })
            print(
                f"{year} {product.upper():9s} {complete:2d}/12 "
                f"ready | missing: {','.join(missing) or '-'} "
                f"| archives awaiting extraction: {','.join(pending) or '-'}"
            )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print("\nDownloader/import manifests (their presence does not prove all requested months succeeded):")
    if not manifests:
        print("  NONE FOUND. Cannot distinguish manually extracted from auto-downloaded files.")
    for p in sorted(manifests):
        with p.open(newline="", encoding="utf-8-sig") as f:
            records = list(csv.DictReader(f))
        src = Counter(x.get("archive_source", "unspecified") for x in records)
        valid = sum(x.get("validated", "").lower() in {"true", "1", "yes"} for x in records)
        print(f"  {p.relative_to(ROOT)}: {len(records)} rows, "
              f"validated={valid}, sources={dict(src)}")
    print(f"\nIgnored unrelated/unrecognised rasters/archives: {ignored_names}")
    print(f"Inventory saved: {OUT}")
    print("Only filenames were checked for availability: run the existing pilot "
          "to validate raster contents/coverage.")


if __name__ == "__main__":
    main()
