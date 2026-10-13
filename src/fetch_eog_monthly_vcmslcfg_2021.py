#!/usr/bin/env python3
"""Download and extract the 2021 monthly EOG VCMSLCFG 75N060W tile.

This downloader uses the current EOG repository layout:
  https://eogdata.mines.edu/nighttime_light/monthly/v10/YYYY/YYYYMM/vcmslcfg/

Authentication
--------------
The script first tests whether direct HTTP access is available. If EOG requires
programmatic authentication, it uses the official OAuth2/OpenID password-grant
flow documented by EOG.

OAuth credentials are never written to disk. They may be supplied as
environment variables:
  EOG_CLIENT_ID
  EOG_CLIENT_SECRET
  EOG_USERNAME
  EOG_PASSWORD

If values are missing, the script prompts interactively. Password and client
secret are entered with getpass and are not echoed.

A short-lived token can alternatively be supplied via EOG_ACCESS_TOKEN, but
OAuth client credentials are preferred because the script can renew tokens
between large downloads.

The script reuses the already-extracted July 2021 pilot files when present,
reuses matching .tgz archives in ~/Downloads, downloads missing months
atomically, extracts the requested layers, validates the raster metadata, and
writes a manifest. Downloaded archives are removed after successful extraction
unless keep_archives is enabled in config.

No credentials are committed to the repository and no annual production data
are modified.
"""

from __future__ import annotations

import argparse
import csv
import getpass
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tarfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen

import rasterio

from project_config import PROJECT_ROOT, load_config


TOKEN_URL = (
    "https://eogauth.mines.edu/realms/eog/protocol/openid-connect/token"
)
ARCHIVE_RE = re.compile(
    r"^SVDNB_npp_(?P<start>\d{8})-(?P<end>\d{8})_"
    r"(?P<tile>[^_]+)_vcmslcfg_v10_.*\.tgz$"
)
USER_AGENT = "lightpollution-research-pipeline/1.0"


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


class EOGAuth:
    def __init__(self):
        self.access_token = os.environ.get("EOG_ACCESS_TOKEN")
        self.client_id = os.environ.get("EOG_CLIENT_ID")
        self.client_secret = os.environ.get("EOG_CLIENT_SECRET")
        self.username = os.environ.get("EOG_USERNAME")
        self.password = os.environ.get("EOG_PASSWORD")
        self.token_obtained_at = 0.0
        self.oauth_ready = all(
            [self.client_id, self.client_secret, self.username, self.password]
        )

    def _prompt_for_oauth(self):
        if self.oauth_ready:
            return
        print(
            "\nEOG direct access is not available. Official OAuth/OpenID "
            "programmatic credentials are required.",
            flush=True,
        )
        print(
            "Credentials are used only in memory for this run and are not "
            "written to disk.",
            flush=True,
        )
        self.client_id = self.client_id or input("EOG client ID: ").strip()
        self.client_secret = self.client_secret or getpass.getpass(
            "EOG client secret: "
        ).strip()
        self.username = self.username or input("EOG username: ").strip()
        self.password = self.password or getpass.getpass(
            "EOG password: "
        )
        self.oauth_ready = all(
            [self.client_id, self.client_secret, self.username, self.password]
        )
        if not self.oauth_ready:
            raise RuntimeError(
                "Incomplete OAuth credentials. EOG browser login is not "
                "automatically shared with command-line downloads."
            )

    def obtain_token(self, *, force: bool = False) -> str:
        if (
            self.access_token
            and not force
            and (not self.oauth_ready or time.monotonic() - self.token_obtained_at < 240)
        ):
            return self.access_token

        if not self.oauth_ready:
            self._prompt_for_oauth()

        data = urlencode(
            {
                "client_id": self.client_id,
                "client_secret": self.client_secret,
                "username": self.username,
                "password": self.password,
                "grant_type": "password",
            }
        ).encode("utf-8")
        req = Request(
            TOKEN_URL,
            data=data,
            method="POST",
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": USER_AGENT,
            },
        )
        try:
            with urlopen(req, timeout=60) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(
                f"EOG token request failed: HTTP {exc.code}: {detail[:500]}"
            ) from exc
        except (URLError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"EOG token request failed: {exc}") from exc

        token = payload.get("access_token")
        if not token:
            raise RuntimeError(
                "EOG token endpoint returned no access_token. "
                f"Response keys: {sorted(payload)}"
            )
        self.access_token = str(token)
        self.token_obtained_at = time.monotonic()
        return self.access_token

    def request(self, url: str, *, bearer: bool, fresh_token: bool = False):
        headers = {"User-Agent": USER_AGENT}
        if bearer:
            token = self.obtain_token(force=fresh_token)
            headers["Authorization"] = f"Bearer {token}"
        return Request(url, headers=headers)


def read_url_text(url: str, auth: EOGAuth, *, bearer: bool) -> str:
    req = auth.request(url, bearer=bearer)
    try:
        with urlopen(req, timeout=60) as response:
            content_type = response.headers.get("Content-Type", "")
            data = response.read()
    except HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code} reading {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"Network error reading {url}: {exc}") from exc

    text = data.decode("utf-8", errors="replace")
    if "text/html" not in content_type.lower() and "<html" not in text.lower():
        raise RuntimeError(
            f"Expected an HTML directory listing at {url}; "
            f"received Content-Type {content_type!r}."
        )
    return text


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
    prefix = f"SVDNB_npp_{year}{month:02d}"
    candidates = []
    for href in parser.hrefs:
        name = Path(href.split("?", 1)[0]).name
        if (
            name.startswith(prefix)
            and f"_{tile}_vcmslcfg_v10_" in name
            and name.endswith(".tgz")
        ):
            candidates.append((name, urljoin(listing_url, href)))

    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected exactly one {tile} VCMSLCFG archive in {listing_url}; "
            f"found {len(candidates)}: {[x[0] for x in candidates]}"
        )
    return candidates[0]


def determine_access_mode(test_url: str, auth: EOGAuth, year: int, month: int, tile: str):
    try:
        html = read_url_text(test_url, auth, bearer=False)
        archive = find_archive_in_listing(
            test_url, html, year=year, month=month, tile=tile
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
            html = read_url_text(test_url, auth, bearer=True)
            archive = find_archive_in_listing(
                test_url, html, year=year, month=month, tile=tile
            )
            print("EOG OAuth/OpenID access: OK", flush=True)
            return True, archive
        except Exception as oauth_exc:
            raise RuntimeError(
                "Unable to access the EOG repository programmatically. "
                "A browser login is a separate session and is not automatically "
                "available to this command-line downloader.\n"
                f"Direct access error: {direct_exc}\n"
                f"OAuth access error: {oauth_exc}"
            ) from oauth_exc


def local_archive_candidates(year: int, month: int, tile: str, archive_dir: Path):
    pattern = f"SVDNB_npp_{year}{month:02d}*_{tile}_vcmslcfg_v10_*.tgz"
    candidates = []
    for root in [archive_dir, Path.home() / "Downloads"]:
        if root.exists():
            candidates.extend(sorted(root.glob(pattern)))
    unique = []
    seen = set()
    for p in candidates:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique.append(p)
    return unique


def existing_extracted_pilot_files(
    year: int,
    month: int,
    tile: str,
    layers: list[str],
):
    pilot = PROJECT_ROOT / "data/raw/viirs/monthly_july_pilot"
    if not pilot.exists():
        return {}
    result = {}
    for layer in layers:
        matches = sorted(
            pilot.glob(
                f"SVDNB_npp_{year}{month:02d}*_{tile}_vcmslcfg_v10_*.{layer}.tif"
            )
        )
        if len(matches) == 1:
            result[layer] = matches[0]
    return result if len(result) == len(layers) else {}


def link_or_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    try:
        os.link(src, dst)
        print(f"  hard-linked existing file: {dst.name}", flush=True)
    except OSError:
        shutil.copy2(src, dst)
        print(f"  copied existing file: {dst.name}", flush=True)


def download_archive(url: str, dst: Path, auth: EOGAuth, *, bearer: bool):
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_suffix(dst.suffix + ".part")
    if part.exists():
        part.unlink()

    req = auth.request(url, bearer=bearer, fresh_token=bearer and auth.oauth_ready)
    print(f"  downloading: {url}", flush=True)
    try:
        with urlopen(req, timeout=120) as response, part.open("wb") as handle:
            total_header = response.headers.get("Content-Length")
            total = int(total_header) if total_header and total_header.isdigit() else None
            written = 0
            last_report = 0
            while True:
                chunk = response.read(8 * 1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                written += len(chunk)
                if written - last_report >= 128 * 1024 * 1024:
                    if total:
                        print(
                            f"    {written / 1024**2:.0f} / "
                            f"{total / 1024**2:.0f} MiB",
                            flush=True,
                        )
                    else:
                        print(
                            f"    {written / 1024**2:.0f} MiB",
                            flush=True,
                        )
                    last_report = written
    except (HTTPError, URLError) as exc:
        if part.exists():
            part.unlink()
        raise RuntimeError(f"Failed downloading {url}: {exc}") from exc

    if part.stat().st_size < 1024:
        part.unlink()
        raise RuntimeError(f"Downloaded archive is implausibly small: {url}")
    if not tarfile.is_tarfile(part):
        part.unlink()
        raise RuntimeError(
            f"Downloaded object is not a valid tar archive: {url}. "
            "This commonly indicates an authentication/login response."
        )

    part.replace(dst)
    print(f"  downloaded {dst.name} ({dst.stat().st_size / 1024**2:.1f} MiB)", flush=True)


def extract_layers(archive: Path, tif_dir: Path, layers: list[str], tile: str):
    tif_dir.mkdir(parents=True, exist_ok=True)
    extracted = {}
    with tarfile.open(archive, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        for layer in layers:
            matches = [
                m
                for m in members
                if Path(m.name).name.endswith(f".{layer}.tif")
                and f"_{tile}_vcmslcfg_v10_" in Path(m.name).name
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


def validate_extracted(paths: dict[str, Path], expected_crs: str, tile: str):
    ref = None
    metadata = {}
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


def parse_archive_dates(name: str):
    m = ARCHIVE_RE.match(name)
    if not m:
        raise RuntimeError(f"Unexpected archive filename: {name}")
    return m.group("start"), m.group("end")


def write_manifest(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "year",
        "month",
        "start_date",
        "end_date",
        "archive_name",
        "archive_source",
        "avg_rade9h_file",
        "cf_cvg_file",
        "cvg_file",
        "width",
        "height",
        "crs",
        "validated",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["downloads"]["eog_monthly_vcmslcfg_2021"]

    year = int(cfg["year"])
    tile = str(cfg["tile"])
    months = [int(m) for m in cfg["months"]]
    layers = [str(x) for x in cfg["layers"]]
    base_url = str(cfg["base_url"]).rstrip("/")
    out_root = PROJECT_ROOT / cfg["output_dir"]
    archive_dir = out_root / "archives"
    tif_dir = out_root / "tif"
    manifest_path = out_root / "manifest.csv"
    keep_archives = bool(cfg.get("keep_archives", False))

    required_layers = {"avg_rade9h", "cf_cvg", "cvg"}
    if set(layers) != required_layers:
        raise ValueError(
            "This downloader currently expects layers "
            f"{sorted(required_layers)}; configured {sorted(layers)}."
        )

    out_root.mkdir(parents=True, exist_ok=True)
    archive_dir.mkdir(parents=True, exist_ok=True)
    tif_dir.mkdir(parents=True, exist_ok=True)

    # Reuse the existing July pilot without consuming additional disk space.
    july_files = existing_extracted_pilot_files(year, 7, tile, layers)
    if july_files:
        for layer, src in july_files.items():
            link_or_copy(src, tif_dir / src.name)

    auth = EOGAuth()
    first_month = months[0]
    first_listing = (
        f"{base_url}/{year}/{year}{first_month:02d}/vcmslcfg/"
    )
    bearer, first_archive = determine_access_mode(
        first_listing,
        auth,
        year,
        first_month,
        tile,
    )

    rows = []
    for month in months:
        print(f"\n===== {year}-{month:02d} VCMSLCFG {tile} =====", flush=True)

        # If all requested layers already exist canonically, validate and record
        # without touching the network for this month.
        existing = {}
        for layer in layers:
            matches = sorted(
                tif_dir.glob(
                    f"SVDNB_npp_{year}{month:02d}*_{tile}_vcmslcfg_v10_*.{layer}.tif"
                )
            )
            if len(matches) == 1:
                existing[layer] = matches[0]
        if len(existing) == len(layers):
            meta = validate_extracted(
                existing,
                config["nuts"]["crs"],
                tile,
            )
            archive_name = re.sub(
                r"\.avg_rade9h\.tif$",
                ".tgz",
                existing["avg_rade9h"].name,
            )
            start, end = parse_archive_dates(archive_name)
            rows.append(
                {
                    "year": year,
                    "month": month,
                    "start_date": start,
                    "end_date": end,
                    "archive_name": archive_name,
                    "archive_source": "existing_extracted",
                    "avg_rade9h_file": existing["avg_rade9h"].name,
                    "cf_cvg_file": existing["cf_cvg"].name,
                    "cvg_file": existing["cvg"].name,
                    "width": meta["avg_rade9h"]["width"],
                    "height": meta["avg_rade9h"]["height"],
                    "crs": meta["avg_rade9h"]["crs"],
                    "validated": True,
                }
            )
            write_manifest(manifest_path, rows)
            print("  already extracted and validated; skipping download", flush=True)
            continue

        local_candidates = local_archive_candidates(
            year, month, tile, archive_dir
        )
        if len(local_candidates) > 1:
            raise RuntimeError(
                f"More than one local archive candidate for {year}-{month:02d}: "
                f"{[str(x) for x in local_candidates]}"
            )

        if local_candidates:
            archive = local_candidates[0]
            archive_source = "local_existing"
            print(f"  reusing local archive: {archive}", flush=True)
        else:
            listing_url = f"{base_url}/{year}/{year}{month:02d}/vcmslcfg/"
            if month == first_month:
                archive_name, archive_url = first_archive
            else:
                html = read_url_text(listing_url, auth, bearer=bearer)
                archive_name, archive_url = find_archive_in_listing(
                    listing_url,
                    html,
                    year=year,
                    month=month,
                    tile=tile,
                )
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
        meta = validate_extracted(
            extracted,
            config["nuts"]["crs"],
            tile,
        )

        start, end = parse_archive_dates(archive.name)
        rows.append(
            {
                "year": year,
                "month": month,
                "start_date": start,
                "end_date": end,
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
        )
        write_manifest(manifest_path, rows)

        if (
            archive_source == "downloaded"
            and not keep_archives
            and archive.parent.resolve() == archive_dir.resolve()
        ):
            archive.unlink()
            print("  removed downloaded archive after validated extraction", flush=True)

    if len(rows) != len(months):
        raise RuntimeError(
            f"Expected {len(months)} completed months; manifest has {len(rows)}."
        )

    print("\nEOG 2021 MONTHLY VCMSLCFG DOWNLOAD COMPLETE")
    print(f"Months: {len(rows)}")
    print(f"Extracted layers: {', '.join(layers)}")
    print(f"Raster directory: {tif_dir}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
