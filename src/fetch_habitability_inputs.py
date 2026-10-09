#!/usr/bin/env python3
"""Download static terrain inputs for the European NUTS3 habitability metric H.

Inputs:
- Copernicus DEM GLO-90 Cloud Optimized GeoTIFF tiles from the public AWS
  bucket (no account required).
- Copernicus Global Land Cover 100 m Collection 3 (2015 base epoch)
  permanent-water and permanent-snow/ice cover-fraction layers from Zenodo.

Only DEM tiles that intersect at least one configured NUTS3 polygon are
downloaded.  wget is used when available, matching the rest of this project.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import urllib.request
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

from project_config import configured_path, load_config


TILE_RE = re.compile(
    r"(Copernicus_DSM_COG_30_([NS])(\d{2})_00_([EW])(\d{3})_00_DEM)"
)


def download(url: str, destination: Path, force: bool, allow_failure: bool = False) -> bool:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not force:
        print(f"Already present: {destination}")
        return True

    if destination.exists():
        destination.unlink()

    try:
        if shutil.which("wget"):
            command = [
                "wget",
                "--continue",
                "--tries=4",
                "--timeout=60",
                "--output-document",
                str(destination),
                url,
            ]
            print("+", " ".join(command), flush=True)
            subprocess.run(command, check=True)
        else:
            print(f"wget not found; urllib download: {url}", flush=True)
            with urllib.request.urlopen(url) as response, destination.open("wb") as out:
                shutil.copyfileobj(response, out)
    except Exception:
        if destination.exists():
            destination.unlink()
        if allow_failure:
            return False
        raise
    return True


def nuts_path(config: dict, release: int) -> Path:
    name = config["nuts"]["filename_template"].format(release=release)
    return configured_path(config, "nuts_raw") / name


def load_all_nuts(config: dict) -> gpd.GeoDataFrame:
    frames = []
    for release in config["nuts"]["releases"]:
        path = nuts_path(config, int(release))
        if not path.exists():
            raise FileNotFoundError(
                f"Missing NUTS {release} boundaries: {path}. "
                "Run: python src/run_pipeline.py nuts"
            )
        g = gpd.read_file(path)[["NUTS_ID", "geometry"]].copy()
        if g.crs is None:
            raise ValueError(f"NUTS file has no CRS: {path}")
        g = g.to_crs("EPSG:4326")
        g["nuts_release"] = int(release)
        frames.append(g)
    return gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True),
        geometry="geometry",
        crs="EPSG:4326",
    )


def parse_tile_name(text: str):
    m = TILE_RE.search(text)
    if not m:
        return None
    name, ns, lat_text, ew, lon_text = m.groups()
    lat = int(lat_text) * (1 if ns == "N" else -1)
    lon = int(lon_text) * (1 if ew == "E" else -1)
    return name, lat, lon


def select_intersecting_tiles(tile_list: Path, nuts: gpd.GeoDataFrame) -> pd.DataFrame:
    unique = {}
    with tile_list.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parsed = parse_tile_name(line)
            if parsed is None:
                continue
            name, lat, lon = parsed
            unique[name] = (lat, lon)

    if not unique:
        raise RuntimeError(
            f"No Copernicus DEM tile names could be parsed from {tile_list}"
        )

    sindex = nuts.sindex
    rows = []
    for name, (lat, lon) in unique.items():
        footprint = box(lon, lat, lon + 1.0, lat + 1.0)
        candidates = sindex.query(footprint, predicate="intersects")
        if len(candidates) == 0:
            continue
        rows.append(
            {
                "tile_id": name,
                "south": lat,
                "west": lon,
                "north": lat + 1.0,
                "east": lon + 1.0,
                "n_intersecting_nuts_features": int(len(candidates)),
            }
        )

    return pd.DataFrame(rows).sort_values(["south", "west"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["habitability"]
    root = configured_path(config, "habitability_raw")
    dem_dir = root / "copernicus_dem_glo90"
    lc_dir = root / "copernicus_landcover_2015"
    root.mkdir(parents=True, exist_ok=True)

    nuts = load_all_nuts(config)

    tile_list_url = cfg["dem"]["tile_list_url"]
    tile_list_path = dem_dir / "tileList.txt"
    download(tile_list_url, tile_list_path, args.force_download)

    manifest = select_intersecting_tiles(tile_list_path, nuts)
    print(
        f"Copernicus GLO-90 tiles intersecting configured NUTS3 releases: "
        f"{len(manifest):,}",
        flush=True,
    )

    base = cfg["dem"]["http_base"].rstrip("/")
    tile_paths = []
    for i, row in manifest.iterrows():
        tile_id = row["tile_id"]
        destination = dem_dir / "tiles" / f"{tile_id}.tif"
        url = f"{base}/{tile_id}/{tile_id}.tif"
        print(f"[{i+1:,}/{len(manifest):,}] {tile_id}", flush=True)
        ok = download(url, destination, args.force_download, allow_failure=False)
        if ok:
            tile_paths.append(str(destination.relative_to(root)))

    manifest["relative_path"] = [
        str((Path("copernicus_dem_glo90") / "tiles" / f"{x}.tif"))
        for x in manifest["tile_id"]
    ]
    manifest["source_url"] = [
        f"{base}/{x}/{x}.tif" for x in manifest["tile_id"]
    ]
    manifest_path = root / "copernicus_dem_glo90_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    for key in ["permanent_water", "snow_ice"]:
        source = cfg["landcover"][key]
        destination = lc_dir / source["filename"]
        download(source["url"], destination, args.force_download)

    print()
    print(f"DEM tile manifest: {manifest_path}")
    print(f"DEM tiles downloaded: {len(tile_paths):,}")
    print(
        "Permanent-water layer: "
        f"{lc_dir / cfg['landcover']['permanent_water']['filename']}"
    )
    print(
        "Snow/ice layer: "
        f"{lc_dir / cfg['landcover']['snow_ice']['filename']}"
    )


if __name__ == "__main__":
    main()
