#!/usr/bin/env python3
"""Download the NUTS3 boundary releases needed by the configured analysis.

Uses wget when available and retains the official GISCO GeoJSON files under
data/raw/nuts. Existing files are left untouched unless --force-download is
supplied.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import urllib.request
from pathlib import Path

from project_config import configured_path, load_config


GISCO_BASE = "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson"


def download(url: str, destination: Path, force: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and not force:
        print(f"Already present: {destination.name}")
        return

    if shutil.which("wget"):
        command = [
            "wget",
            "--continue",
            "--output-document",
            str(destination),
            url,
        ]
        print("+", " ".join(command), flush=True)
        subprocess.run(command, check=True)
    else:
        print("wget not found; falling back to Python urllib.", flush=True)
        with urllib.request.urlopen(url) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    root = configured_path(config, "nuts_raw")
    template = config["nuts"]["filename_template"]

    for release in config["nuts"]["releases"]:
        filename = template.format(release=release)
        url = f"{GISCO_BASE}/{filename}"
        download(url, root / filename, args.force_download)

    print()
    print("Required NUTS releases are present:")
    for release in config["nuts"]["releases"]:
        filename = template.format(release=release)
        path = root / filename
        print(f"  {release}: {path}")


if __name__ == "__main__":
    main()
