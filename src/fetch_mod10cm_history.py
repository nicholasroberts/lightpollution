#!/usr/bin/env python3
"""Discover and download historical MOD10CM Collection 6.1 monthly snow files.

The downloader uses NASA CMR for discovery and Earthdata Login for protected
file access.  It is deliberately idempotent: already-downloaded HDF files are
reused, partial downloads are resumed, and an inventory manifest is written
after every successful run.

Authentication is never stored in the repository.  Two modes are supported:

1. Preferred for unattended runs:
       export EARTHDATA_TOKEN='...'
   where the value is an Earthdata Login *user token* (not an AppEEARS API
   session token).

2. Standard Earthdata .netrc authentication:
       ~/.netrc
   containing a machine entry for urs.earthdata.nasa.gov.

The stage downloads all configured months so the historical snow archive can be
reused for later seasonal diagnostics.  The subsequent historical surface-state
model can still pre-specify January-February as the primary control window.
"""

from __future__ import annotations

import argparse
import calendar
import csv
import hashlib
import json
import netrc
import os
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

try:
    import earthaccess
except ImportError as exc:
    raise RuntimeError(
        "The historical MOD10CM downloader now uses NASA's earthaccess "
        "library for authenticated Earthdata Cloud downloads. Install the "
        "project requirements with: python -m pip install -r requirements.txt"
    ) from exc

from project_config import PROJECT_ROOT, load_config


CMR_URL = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
HDF_RE = re.compile(
    r"^MOD10CM\.A(?P<year>\d{4})(?P<doy>\d{3})\.061\..*\.hdf$"
)


def month_from_year_doy(year: int, doy: int) -> int:
    import datetime as dt

    return (dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)).month


def cmr_search(cfg: dict) -> list[dict[str, Any]]:
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    params = {
        "short_name": str(cfg["short_name"]),
        "version": str(cfg["version"]),
        "temporal": (
            f"{first_year}-01-01T00:00:00Z,"
            f"{last_year}-12-31T23:59:59Z"
        ),
        "downloadable": "true",
        "page_size": "2000",
    }
    url = CMR_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "lightpollution-MOD10CM-downloader/1.0",
            "Accept": "application/vnd.nasa.cmr.umm_results+json",
        },
    )
    print(f"Querying CMR: {url}", flush=True)
    with urllib.request.urlopen(req, timeout=120) as response:
        payload = json.load(response)

    items = payload.get("items", [])
    hits = int(payload.get("hits", len(items)))
    if hits > len(items):
        raise RuntimeError(
            f"CMR returned {len(items)} of {hits} hits. Increase page_size or "
            "add pagination before continuing."
        )
    if not items:
        raise RuntimeError("CMR returned no MOD10CM granules.")
    return items


def checksum_from_umm(umm: dict[str, Any]) -> tuple[str | None, str | None]:
    dg = umm.get("DataGranule") or {}
    checksum = dg.get("Checksum") or {}
    value = checksum.get("Value")
    algorithm = checksum.get("Algorithm")
    if value and algorithm:
        return str(algorithm).upper(), str(value).lower()
    return None, None


def choose_data_url(umm: dict[str, Any]) -> str:
    """Choose an HTTP(S) HDF download URL from a CMR granule record.

    CMR can advertise both Earthdata Cloud S3 object URLs and ordinary HTTPS
    distribution URLs for the same granule.  curl cannot read s3:// URLs unless
    built with S3 support/extra tooling, and the repository deliberately uses
    standard curl + Earthdata authentication.  Therefore S3 URLs are *not*
    candidates here even when CMR labels them "GET DATA VIA DIRECT ACCESS".
    """
    candidates: list[tuple[int, str]] = []
    seen_related: list[str] = []

    for rel in umm.get("RelatedUrls") or []:
        url = str(rel.get("URL") or "")
        if not url:
            continue
        seen_related.append(url)

        parsed = urllib.parse.urlparse(url)
        if parsed.scheme not in ("http", "https"):
            continue

        name = Path(parsed.path).name
        if not HDF_RE.match(name):
            continue

        rtype = str(rel.get("Type") or "").upper()
        subtype = str(rel.get("Subtype") or "").upper()

        # We require an ordinary web URL.  Within that constraint prefer
        # explicit GET DATA links and HTTPS.
        score = 0
        if rtype == "GET DATA":
            score += 30
        if "DIRECT" in subtype:
            score += 5
        if parsed.scheme == "https":
            score += 10
        candidates.append((score, url))

    if not candidates:
        preview = "\n".join(f"  {u}" for u in seen_related[:12])
        raise RuntimeError(
            f"No HTTP(S) HDF data URL found in CMR record for "
            f"{umm.get('GranuleUR', '<unknown granule>')}. "
            "CMR related URLs included:\n"
            f"{preview}"
        )

    candidates.sort(key=lambda x: (-x[0], x[1]))
    return candidates[0][1]

def build_inventory(items: list[dict[str, Any]], cfg: dict) -> list[dict[str, Any]]:
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    requested_months = {int(x) for x in cfg["months"]}
    by_key: dict[tuple[int, int], dict[str, Any]] = {}

    for item in items:
        umm = item.get("umm") or {}
        try:
            url = choose_data_url(umm)
        except RuntimeError:
            continue
        name = Path(urllib.parse.urlparse(url).path).name
        match = HDF_RE.match(name)
        if not match:
            continue
        year = int(match.group("year"))
        month = month_from_year_doy(year, int(match.group("doy")))
        if not (first_year <= year <= last_year) or month not in requested_months:
            continue
        key = (year, month)
        if key in by_key:
            raise RuntimeError(
                f"CMR returned duplicate MOD10CM granules for {year}-{month:02d}: "
                f"{by_key[key]['filename']} and {name}"
            )
        algorithm, checksum = checksum_from_umm(umm)
        by_key[key] = {
            "year": year,
            "month": month,
            "filename": name,
            "url": url,
            "checksum_algorithm": algorithm or "",
            "checksum_expected": checksum or "",
        }

    expected = {
        (year, month)
        for year in range(first_year, last_year + 1)
        for month in requested_months
    }
    allowed_missing = {
        tuple(int(v) for v in str(x).split("-"))
        for x in cfg.get("allowed_missing_year_month", [])
    }
    missing = sorted(expected - set(by_key))
    unexpected_missing = [x for x in missing if x not in allowed_missing]
    if unexpected_missing:
        raise RuntimeError(
            "CMR inventory is missing requested MOD10CM months not listed as "
            f"allowed missing: {unexpected_missing}"
        )
    if missing:
        print(
            "CMR missing allowed month(s): "
            + ", ".join(f"{y}-{m:02d}" for y, m in missing),
            flush=True,
        )

    return [by_key[key] for key in sorted(by_key)]


def authenticate_earthaccess():
    """Authenticate through NASA's supported earthaccess client.

    earthaccess handles the OAuth/EDL session needed by Earthdata Cloud HTTPS
    endpoints and is more robust than reproducing the redirect/cookie flow with
    raw curl.  The user's ~/.netrc remains the credential source.
    """
    netrc_path = Path.home() / ".netrc"
    if not netrc_path.exists():
        raise RuntimeError(
            "No ~/.netrc found. Configure Earthdata Login credentials first."
        )
    mode = netrc_path.stat().st_mode & 0o777
    if mode & 0o077:
        raise RuntimeError(
            f"{netrc_path} permissions are {oct(mode)}; run: chmod 600 ~/.netrc"
        )

    try:
        auth = earthaccess.login(strategy="netrc")
    except Exception as exc:
        raise RuntimeError(
            "earthaccess could not authenticate with Earthdata Login using "
            "~/.netrc. Check the EDL username/password and any required "
            "application authorization."
        ) from exc

    if not getattr(auth, "authenticated", False):
        raise RuntimeError(
            "earthaccess returned without an authenticated Earthdata session."
        )
    return auth

def compute_checksum(path: Path, algorithm: str) -> str:
    name = algorithm.replace("-", "").lower()
    supported = {"md5": "md5", "sha1": "sha1", "sha256": "sha256"}
    if name not in supported:
        return ""
    h = hashlib.new(supported[name])
    with path.open("rb") as src:
        for chunk in iter(lambda: src.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest().lower()


def validate_existing(
    path: Path, checksum_algorithm: str, checksum_expected: str
) -> tuple[bool, str]:
    if not path.exists() or path.stat().st_size <= 1024:
        return False, ""
    if checksum_algorithm and checksum_expected:
        actual = compute_checksum(path, checksum_algorithm)
        if actual:
            return actual == checksum_expected.lower(), actual
    return True, ""


def download_one(
    row: dict[str, Any],
    raw_dir: Path,
) -> dict[str, Any]:
    target = raw_dir / row["filename"]
    ok, checksum_actual = validate_existing(
        target,
        row["checksum_algorithm"],
        row["checksum_expected"],
    )
    if ok:
        status = "reused"
        print(f"REUSE {row['year']}-{row['month']:02d} {target.name}", flush=True)
    else:
        if target.exists():
            print(f"Existing file failed validation; redownloading {target.name}")
            target.unlink()

        print(
            f"DOWNLOAD {row['year']}-{row['month']:02d} {target.name}",
            flush=True,
        )

        # earthaccess accepts HTTP URL strings when the provider is supplied.
        # Running one granule at a time preserves the existing per-file
        # checksum, progress-manifest, and clean-restart semantics.
        try:
            paths = earthaccess.download(
                [row["url"]],
                local_path=str(raw_dir),
                provider="NSIDC_CPRD",
                threads=1,
                show_progress=True,
                force=True,
            )
        except Exception as exc:
            raise RuntimeError(
                f"earthaccess download failed for "
                f"{row['year']}-{row['month']:02d}: {row['url']}"
            ) from exc

        if not paths:
            raise RuntimeError(
                f"earthaccess returned no local file for {row['filename']}"
            )

        # earthaccess normally preserves the source filename.  Resolve its
        # returned path explicitly so a future library change cannot silently
        # leave the expected archive entry absent.
        returned = Path(paths[0])
        if not returned.is_absolute():
            returned = (PROJECT_ROOT / returned).resolve()
        if not target.exists() and returned.exists():
            if returned.name != target.name:
                raise RuntimeError(
                    f"earthaccess downloaded {returned.name}, expected "
                    f"{target.name}"
                )
            if returned.resolve() != target.resolve():
                returned.replace(target)

        if not target.exists() or target.stat().st_size <= 1024:
            raise RuntimeError(
                f"Downloaded file is missing or unexpectedly small: {target}"
            )

        ok, checksum_actual = validate_existing(
            target,
            row["checksum_algorithm"],
            row["checksum_expected"],
        )
        if not ok:
            raise RuntimeError(
                f"Checksum validation failed after download: {target.name}"
            )
        status = "downloaded"

    out = dict(row)
    out.update(
        {
            "local_path": str(target.relative_to(PROJECT_ROOT)),
            "size_bytes": int(target.stat().st_size),
            "checksum_actual": checksum_actual,
            "status": status,
        }
    )
    return out

def write_manifest(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "year", "month", "filename", "url", "local_path", "size_bytes",
        "checksum_algorithm", "checksum_expected", "checksum_actual", "status",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def inventory_summary(
    rows: list[dict[str, Any]], cfg: dict, manifest: Path
) -> None:
    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    months = sorted(int(x) for x in cfg["months"])
    expected_total = (last_year - first_year + 1) * len(months)
    allowed_missing = len(cfg.get("allowed_missing_year_month", []))
    downloaded = sum(r["status"] == "downloaded" for r in rows)
    reused = sum(r["status"] == "reused" for r in rows)
    size_gib = sum(int(r["size_bytes"]) for r in rows) / 1024**3

    print("\nMOD10CM HISTORICAL DOWNLOAD COMPLETE")
    print(f"Years: {first_year}-{last_year}")
    print(f"Months requested per year: {months}")
    print(
        f"Available granules: {len(rows)} / nominal {expected_total} "
        f"(configured allowed-missing entries: {allowed_missing})"
    )
    print(f"Downloaded this run: {downloaded}")
    print(f"Reused existing: {reused}")
    print(f"Local archive size: {size_gib:.2f} GiB")
    print(f"Manifest: {manifest}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["mod10cm_history_download"]
    raw_dir = PROJECT_ROOT / str(cfg["raw_dir"])
    raw_dir.mkdir(parents=True, exist_ok=True)

    items = cmr_search(cfg)
    inventory = build_inventory(items, cfg)
    print(f"CMR discovered {len(inventory)} requested MOD10CM granules.", flush=True)

    authenticate_earthaccess()
    print("Earthdata authentication mode: earthaccess/netrc", flush=True)

    completed: list[dict[str, Any]] = []
    manifest = PROJECT_ROOT / str(cfg["manifest_file"])
    try:
        for row in inventory:
            completed.append(
                download_one(row, raw_dir)
            )
            # Persist progress so an interrupted long run remains auditable.
            write_manifest(completed, manifest)
    except KeyboardInterrupt:
        write_manifest(completed, manifest)
        print("\nInterrupted; completed manifest preserved.", file=sys.stderr)
        raise

    write_manifest(completed, manifest)
    inventory_summary(completed, cfg, manifest)


if __name__ == "__main__":
    main()
