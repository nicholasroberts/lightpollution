#!/usr/bin/env python3
"""Map and diagnose the 2021 VIIRS anomaly using monthly VCMCFG/VCMSLCFG data.

This is a diagnostic pipeline stage, not a replacement for the annual production
series. It asks whether the 2021 annual radiance pulse is concentrated in a few
months and/or locations, and whether it behaves differently in the monthly
VCMCFG (stray-light observations excluded) and VCMSLCFG (stray-light corrected)
products.

Expected source data
--------------------
The script scans the existing VIIRS raw-data tree recursively:

    data/raw/viirs/

This accommodates the current project layout, including directories such as
monthly_vcmslcfg_2020, monthly_vcmslcfg_2021 and monthly_vcmslcfg_2022, while
leaving the annual global VIIRS files in the same raw-data root. For the full
comparison it needs all months of 2020, 2021 and 2022 for both vcmcfg and
vcmslcfg.
For each product-month the script needs:

* the average-radiance GeoTIFF (filename containing avg_rade9); and
* the cloud-free coverage GeoTIFF (filename containing cf_cvg).

The current EOG site requires authenticated access. Downloading is intentionally
kept outside this script so credentials are never stored in the repository.
Non-tiled/global monthly files are preferred because the locked NUTS2016 cohort
can include overseas NUTS regions. Tiled inputs are accepted only when exactly
one candidate covers the requested product-month; the script will not silently
mosaic incomplete tile sets.

Locked analytical design
------------------------
* fixed NUTS2016 geography;
* exact 1,514-region cohort from the existing fixed-NUTS2016 analysis;
* fixed population-density classes at 30/100/300/1000 people km^-2;
* 2021 compared month-by-month with mean(2020, 2022) for the same calendar month;
* cloud-free coverage (cf_cvg > 0) is used as a validity mask;
* NUTS3 monthly means use exact polygon overlap and spherical physical-area
  weighting, matching the production zonal-statistics design;
* region-month summaries used for inference require at least 50% valid NUTS3
  area by default (configurable with --min-valid-area-fraction).

Outputs
-------
Processed tables are written under:

    data/processed/monthly_2021_anomaly/

and figures under:

    figures/annual_change/monthly_2021_anomaly/

The extraction is checkpointed by product-month so interrupted runs can resume.
"""

from __future__ import annotations

import argparse
import calendar
import gc
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterSource
from matplotlib.colors import SymLogNorm, TwoSlopeNorm
from rasterio.enums import Resampling
from rasterio.transform import from_bounds
from rasterio.warp import reproject
from rasterio.windows import Window
from shapely.geometry import box

from diagnose_fixed_nuts2016_2021_reversal import load_locked_panel
from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


YEARS = (2020, 2021, 2022)
PRODUCTS = ("vcmcfg", "vcmslcfg")
DENSITY_ORDER = ("<30", "30-100", "100-300", "300-1000", ">=1000")
EXPECTED_COHORT_SIZE = 1514
DEFAULT_MIN_VALID_AREA_FRACTION = 0.50
DEFAULT_MAP_BOUNDS = (-25.0, 34.0, 45.0, 72.5)
DEFAULT_MAX_PIXELS_ACROSS = 1400
DEFAULT_PERCENT_MAP_MIN_RADIANCE = 0.20

PRODUCT_LABELS = {
    "vcmcfg": "VCMCFG (stray-light observations excluded)",
    "vcmslcfg": "VCMSLCFG (stray-light corrected)",
}

YM_RE = re.compile(r"(?<!\d)(20\d{2})(0[1-9]|1[0-2])\d{2}(?!\d)")


@dataclass(frozen=True)
class RasterPair:
    product: str
    year: int
    month: int
    radiance: Path
    cloudfree: Path

    @property
    def yyyymm(self) -> str:
        return f"{self.year}{self.month:02d}"


class CoverageMaskedRasterSource(RasterSource):
    """RasterSource that masks avg radiance wherever cf_cvg <= 0."""

    def __init__(self, radiance_path: Path, cloudfree_path: Path, *, name: str):
        super().__init__()
        self.radiance_path = Path(radiance_path)
        self.cloudfree_path = Path(cloudfree_path)
        self.rad = rasterio.open(self.radiance_path)
        self.cov = rasterio.open(self.cloudfree_path)
        self._validate_grids()
        self.set_name(name)

    def _validate_grids(self) -> None:
        if self.rad.crs != self.cov.crs:
            raise ValueError(
                f"CRS mismatch:\n  {self.radiance_path}\n  {self.cloudfree_path}"
            )
        if self.rad.transform != self.cov.transform:
            raise ValueError(
                f"Grid transform mismatch:\n  {self.radiance_path}\n  {self.cloudfree_path}"
            )
        if (self.rad.width, self.rad.height) != (self.cov.width, self.cov.height):
            raise ValueError(
                f"Raster shape mismatch:\n  {self.radiance_path}\n  {self.cloudfree_path}"
            )
        if self.rad.count != 1 or self.cov.count != 1:
            raise ValueError("Monthly VIIRS inputs must each contain exactly one band.")
        if self.rad.crs is None:
            raise ValueError(f"Missing CRS: {self.radiance_path}")

    def res(self):
        return (abs(self.rad.transform.a), abs(self.rad.transform.e))

    def extent(self):
        b = self.rad.bounds
        return (b.left, b.bottom, b.right, b.top)

    def srs_wkt(self):
        return self.rad.crs.wkt

    def nodata_value(self):
        return None

    def read_window(self, x0, y0, nx, ny):
        if nx == 0 or ny == 0:
            return np.ma.array(np.empty((0, 0), dtype=np.float32), mask=True)
        window = Window(int(x0), int(y0), int(nx), int(ny))
        rad = self.rad.read(1, window=window, masked=True)
        cov = self.cov.read(1, window=window, masked=True)
        rad_data = np.asarray(rad.data, dtype=np.float32)
        cov_data = np.asarray(cov.data)
        mask = (
            np.ma.getmaskarray(rad)
            | np.ma.getmaskarray(cov)
            | ~np.isfinite(rad_data)
            | (cov_data <= 0)
        )
        return np.ma.array(rad_data, mask=mask, copy=False)

    def close(self) -> None:
        self.rad.close()
        self.cov.close()


class ConstantGridRasterSource(RasterSource):
    """Constant-one raster on a reference grid, used for exact polygon area."""

    def __init__(self, reference_path: Path, *, name: str = "constant"):
        super().__init__()
        self.path = Path(reference_path)
        self.ds = rasterio.open(self.path)
        if self.ds.crs is None:
            raise ValueError(f"Missing CRS: {self.path}")
        self.set_name(name)

    def res(self):
        return (abs(self.ds.transform.a), abs(self.ds.transform.e))

    def extent(self):
        b = self.ds.bounds
        return (b.left, b.bottom, b.right, b.top)

    def srs_wkt(self):
        return self.ds.crs.wkt

    def nodata_value(self):
        return None

    def read_window(self, x0, y0, nx, ny):
        return np.ones((int(ny), int(nx)), dtype=np.uint8)

    def close(self) -> None:
        self.ds.close()


def month_name(month: int) -> str:
    return calendar.month_abbr[int(month)]


def _extract_year_month(name: str) -> tuple[int, int] | None:
    """Extract YYYY-MM from EOG paths containing YYYYMMDD date stamps.

    Canonical monthly files are named like SVDNB_npp_20200101-20200131_...
    so the year-month token is the first six digits of an eight-digit date.
    Do not require a non-digit immediately after MM.
    """
    match = YM_RE.search(name)
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def _classify_file(path: Path) -> tuple[str, int, int, str] | None:
    # Use the full path, not only the filename.  The project stores monthly
    # products in product/year-labelled directories, so parent-directory names
    # are useful metadata and allow the inventory to scan data/raw/viirs
    # directly without reorganising existing raw data.
    name = str(path).lower()
    ym = _extract_year_month(name)
    if ym is None:
        return None

    if "vcmslcfg" in name:
        product = "vcmslcfg"
    elif "vcmcfg" in name:
        product = "vcmcfg"
    else:
        return None

    if "cf_cvg" in name:
        layer = "cloudfree"
    elif "avg_rade9" in name or "avg_rad" in name:
        layer = "radiance"
    else:
        return None

    year, month = ym
    return product, year, month, layer


def _candidate_score(path: Path) -> tuple[int, int, str]:
    """Prefer canonical monthly directories and avoid duplicate July pilot files."""
    full = str(path).lower()
    name = path.name.lower()
    if "monthly_july_pilot" in full:
        tier = 0
    elif "global" in name:
        tier = 4
    elif re.search(r"monthly_vcm(?:sl)?cfg_20\\d{2}/tif", full):
        tier = 3
    elif "75n060w" in name:
        tier = 2
    else:
        tier = 1
    return (tier, -len(str(path)), str(path))


def _choose_one(paths: list[Path], *, key: str) -> Path:
    if not paths:
        raise FileNotFoundError(key)
    ranked = sorted(paths, key=_candidate_score, reverse=True)
    top_tier = _candidate_score(ranked[0])[0]
    top = [p for p in ranked if _candidate_score(p)[0] == top_tier]
    if len(top) > 1:
        lines = "\n".join(f"  {p}" for p in top)
        raise RuntimeError(
            f"Ambiguous monthly input for {key}. Multiple equally preferred files:\n"
            f"{lines}\nRemove duplicates or use a dedicated --input-dir."
        )
    return ranked[0]


def discover_inputs(input_dir: Path) -> tuple[dict[tuple[str, int, int], RasterPair], list[str]]:
    buckets: dict[tuple[str, int, int, str], list[Path]] = {}
    if input_dir.exists():
        for path in sorted(input_dir.rglob("*")):
            if not path.is_file() or path.suffix.lower() not in {".tif", ".tiff"}:
                continue
            info = _classify_file(path)
            if info is None:
                continue
            product, year, month, layer = info
            if product not in PRODUCTS or year not in YEARS:
                continue
            buckets.setdefault((product, year, month, layer), []).append(path)

    pairs: dict[tuple[str, int, int], RasterPair] = {}
    missing: list[str] = []

    for product in PRODUCTS:
        for year in YEARS:
            for month in range(1, 13):
                key = (product, year, month)
                rad_paths = buckets.get((*key, "radiance"), [])
                cov_paths = buckets.get((*key, "cloudfree"), [])
                if not rad_paths:
                    missing.append(
                        f"{product} {year}-{month:02d}: average radiance (avg_rade9*)"
                    )
                    continue
                if not cov_paths:
                    missing.append(
                        f"{product} {year}-{month:02d}: cloud-free coverage (cf_cvg)"
                    )
                    continue
                radiance = _choose_one(
                    rad_paths, key=f"{product} {year}-{month:02d} radiance"
                )
                cloudfree = _choose_one(
                    cov_paths, key=f"{product} {year}-{month:02d} cf_cvg"
                )
                pairs[key] = RasterPair(product, year, month, radiance, cloudfree)

    return pairs, missing


def print_inventory(
    pairs: dict[tuple[str, int, int], RasterPair],
    missing: list[str],
    input_dir: Path,
) -> None:
    print("MONTHLY VIIRS INPUT INVENTORY")
    print("=" * 78)
    print(f"Input root: {input_dir}")
    print(f"Complete product-month pairs: {len(pairs)}/72")
    for product in PRODUCTS:
        n = sum(1 for p, _, _ in pairs if p == product)
        print(f"  {product}: {n}/36")
    if missing:
        print(f"\nMissing layers ({len(missing)}):")
        for item in missing:
            print(f"  - {item}")
    else:
        print("\nAll 2020-2022 VCMCFG/VCMSLCFG monthly inputs are present.")


def _first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def assign_density_band(density: pd.Series) -> pd.Series:
    values = pd.to_numeric(density, errors="coerce")
    return pd.cut(
        values,
        bins=[-np.inf, 30, 100, 300, 1000, np.inf],
        labels=list(DENSITY_ORDER),
        right=False,
        ordered=True,
    ).astype("string")


def load_locked_cohort(processed: Path, cohort_file: Path | None) -> pd.DataFrame:
    if cohort_file is None:
        cohort_file = _first_existing(
            [
                processed / "fixed_nuts2016_region_slopes_2014_2024.csv",
                processed / "fixed_nuts2016_region_slopes_2014_2024.parquet",
            ]
        )
    if cohort_file is None or not cohort_file.exists():
        raise FileNotFoundError(
            "Could not find the locked fixed-NUTS2016 cohort file. Expected one of:\n"
            f"  {processed / 'fixed_nuts2016_region_slopes_2014_2024.csv'}\n"
            f"  {processed / 'fixed_nuts2016_region_slopes_2014_2024.parquet'}\n"
            "Use --cohort-file to point to the existing 1,514-region table. "
            "The script will not reconstruct a different cohort silently."
        )

    if cohort_file.suffix.lower() == ".parquet":
        raw = pd.read_parquet(cohort_file)
    else:
        raw = pd.read_csv(cohort_file)

    id_candidates = ["NUTS_ID", "nuts_id", "nuts3", "region_id"]
    id_col = next((c for c in id_candidates if c in raw.columns), None)
    if id_col is None:
        raise ValueError(
            f"No NUTS identifier column found in {cohort_file}. "
            f"Columns: {list(raw.columns)}"
        )

    band_candidates = [
        "density_band",
        "population_density_band",
        "density_class",
        "pop_density_band",
    ]
    band_col = next((c for c in band_candidates if c in raw.columns), None)

    density_candidates = [
        "reference_population_density_median",
        "population_density_median_2014_2020",
        "density_median_2014_2020",
        "reference_density",
        "population_density",
        "density",
    ]
    density_col = next((c for c in density_candidates if c in raw.columns), None)

    cohort = pd.DataFrame({"NUTS_ID": raw[id_col].astype(str)})
    if band_col is not None:
        cohort["density_band"] = raw[band_col].astype("string")
    elif density_col is not None:
        cohort["reference_population_density_median"] = pd.to_numeric(
            raw[density_col], errors="coerce"
        )
        cohort["density_band"] = assign_density_band(
            cohort["reference_population_density_median"]
        )
    else:
        raise ValueError(
            f"Could not find fixed density-band or reference-density column in "
            f"{cohort_file}. Columns: {list(raw.columns)}"
        )

    if density_col is not None and "reference_population_density_median" not in cohort:
        cohort["reference_population_density_median"] = pd.to_numeric(
            raw[density_col], errors="coerce"
        )

    cohort = cohort.drop_duplicates("NUTS_ID").reset_index(drop=True)
    if len(cohort) != EXPECTED_COHORT_SIZE:
        raise ValueError(
            "Locked cohort mismatch: expected exactly "
            f"{EXPECTED_COHORT_SIZE:,} unique NUTS3 regions, found {len(cohort):,} "
            f"in {cohort_file}. Refusing to change cohort implicitly."
        )

    invalid_bands = sorted(set(cohort["density_band"].dropna()) - set(DENSITY_ORDER))
    if invalid_bands:
        raise ValueError(f"Unexpected density-band labels: {invalid_bands}")
    if cohort["density_band"].isna().any():
        raise ValueError("Some locked-cohort regions have no density band.")

    cohort["density_band"] = pd.Categorical(
        cohort["density_band"], categories=list(DENSITY_ORDER), ordered=True
    )

    print(f"Locked cohort: {len(cohort):,} NUTS3 regions from {cohort_file}")
    print(cohort["density_band"].value_counts(sort=False).to_string())
    return cohort


def load_nuts2016(config: dict, cohort: pd.DataFrame) -> gpd.GeoDataFrame:
    nuts_raw = configured_path(config, "nuts_raw")
    template = config["nuts"]["filename_template"]
    path = nuts_raw / template.format(release=2016)
    if not path.exists():
        raise FileNotFoundError(
            f"NUTS2016 level-3 geometry not found: {path}\n"
            "Run the existing NUTS boundary stage first."
        )

    nuts = gpd.read_file(path)
    if nuts.crs is None:
        raise ValueError(f"Missing CRS: {path}")
    if nuts.crs.to_string() != "EPSG:4326":
        nuts = nuts.to_crs("EPSG:4326")

    id_field = config["nuts"]["id_field"]
    keep_cols = [
        c
        for c in [id_field, "CNTR_CODE", "NUTS_NAME", "NAME_LATN", "geometry"]
        if c in nuts.columns
    ]
    nuts = nuts[keep_cols].copy()
    nuts = nuts.rename(columns={id_field: "NUTS_ID"})
    nuts["NUTS_ID"] = nuts["NUTS_ID"].astype(str)
    nuts = nuts.merge(cohort, on="NUTS_ID", how="inner", validate="one_to_one")

    if len(nuts) != EXPECTED_COHORT_SIZE:
        missing = sorted(set(cohort["NUTS_ID"]) - set(nuts["NUTS_ID"]))
        raise ValueError(
            f"Only {len(nuts):,}/{EXPECTED_COHORT_SIZE:,} locked regions found "
            f"in NUTS2016 geometry. First missing IDs: {missing[:20]}"
        )

    if "CNTR_CODE" not in nuts.columns:
        nuts["CNTR_CODE"] = nuts["NUTS_ID"].str[:2]

    return nuts


def exact_total_region_area(
    nuts: gpd.GeoDataFrame, reference_raster: Path
) -> pd.DataFrame:
    source = ConstantGridRasterSource(reference_raster, name="grid")
    try:
        area = exact_extract(
            source,
            nuts,
            ["region_area_km2=count(coverage_weight=area_spherical_km2)"],
            include_cols=["NUTS_ID"],
            output="pandas",
            strategy="raster-sequential",
            progress=False,
        )
    finally:
        source.close()

    if area["region_area_km2"].isna().any() or (area["region_area_km2"] <= 0).any():
        raise ValueError(
            "Failed to calculate positive exact area for every locked NUTS3 region."
        )
    return area[["NUTS_ID", "region_area_km2"]]


def extract_product_month(
    pair: RasterPair,
    nuts: gpd.GeoDataFrame,
    region_area: pd.DataFrame,
    min_valid_fraction: float,
) -> pd.DataFrame:
    print(
        f"Extracting {pair.product} {pair.year}-{pair.month:02d}: "
        f"{pair.radiance.name}",
        flush=True,
    )

    source = CoverageMaskedRasterSource(
        pair.radiance, pair.cloudfree, name=f"{pair.product}_{pair.yyyymm}"
    )
    try:
        result = exact_extract(
            source,
            nuts,
            [
                "radiance_mean=mean(coverage_weight=area_spherical_m2)",
                "valid_area_km2=count(coverage_weight=area_spherical_km2)",
            ],
            include_cols=["NUTS_ID", "CNTR_CODE", "density_band"],
            output="pandas",
            strategy="raster-sequential",
            progress=True,
        )
    finally:
        source.close()

    result = result.merge(region_area, on="NUTS_ID", how="left", validate="one_to_one")
    result["valid_area_fraction"] = (
        result["valid_area_km2"] / result["region_area_km2"]
    ).clip(lower=0.0, upper=1.0)
    result["analysis_valid"] = (
        np.isfinite(result["radiance_mean"])
        & (result["valid_area_fraction"] >= min_valid_fraction)
    )
    result["product"] = pair.product
    result["year"] = pair.year
    result["month"] = pair.month
    result["yyyymm"] = int(pair.yyyymm)
    result["radiance_file"] = pair.radiance.name
    result["cloudfree_file"] = pair.cloudfree.name
    result["min_valid_area_fraction"] = min_valid_fraction
    result["density_band"] = pd.Categorical(
        result["density_band"], categories=list(DENSITY_ORDER), ordered=True
    )
    return result


def run_zonal_extraction(
    pairs: dict[tuple[str, int, int], RasterPair],
    nuts: gpd.GeoDataFrame,
    processed_out: Path,
    min_valid_fraction: float,
    force: bool,
) -> pd.DataFrame:
    zonal_dir = processed_out / "zonal"
    zonal_dir.mkdir(parents=True, exist_ok=True)

    reference = pairs[(PRODUCTS[0], YEARS[0], 1)].radiance
    region_area_path = processed_out / "locked_nuts2016_exact_area.parquet"
    if region_area_path.exists() and not force:
        region_area = pd.read_parquet(region_area_path)
    else:
        print("Calculating exact locked-cohort NUTS3 area on the monthly VIIRS grid...")
        region_area = exact_total_region_area(nuts, reference)
        region_area.to_parquet(region_area_path, index=False)

    frames: list[pd.DataFrame] = []
    for product in PRODUCTS:
        for year in YEARS:
            for month in range(1, 13):
                pair = pairs[(product, year, month)]
                checkpoint = zonal_dir / f"{product}_{year}{month:02d}.parquet"

                if checkpoint.exists() and not force:
                    print(f"Using checkpoint: {checkpoint}")
                    frame = pd.read_parquet(checkpoint)
                    frame["analysis_valid"] = (
                        np.isfinite(frame["radiance_mean"])
                        & (frame["valid_area_fraction"] >= min_valid_fraction)
                    )
                    frame["min_valid_area_fraction"] = min_valid_fraction
                else:
                    frame = extract_product_month(
                        pair, nuts, region_area, min_valid_fraction
                    )
                    frame.to_parquet(checkpoint, index=False)
                    print(f"Wrote checkpoint: {checkpoint}")

                frames.append(frame)
                gc.collect()

    combined = pd.concat(frames, ignore_index=True)
    expected_rows = EXPECTED_COHORT_SIZE * len(PRODUCTS) * len(YEARS) * 12
    if len(combined) != expected_rows:
        raise ValueError(
            f"Unexpected combined row count: {len(combined):,}; "
            f"expected {expected_rows:,}"
        )

    combined_path = processed_out / "monthly_nuts2016_radiance_2020_2022.parquet"
    combined.to_parquet(combined_path, index=False)
    combined.to_csv(combined_path.with_suffix(".csv"), index=False)
    print(f"Wrote combined monthly zonal table: {combined_path}")
    return combined


def _geometric_pct(ratio: pd.Series) -> float:
    vals = pd.to_numeric(ratio, errors="coerce")
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.empty:
        return np.nan
    return float((np.exp(np.log(vals).mean()) - 1.0) * 100.0)


def build_region_anomalies(monthly: pd.DataFrame) -> pd.DataFrame:
    use = monthly.copy()
    use["radiance_analysis"] = use["radiance_mean"].where(use["analysis_valid"])

    id_cols = ["product", "month", "NUTS_ID", "CNTR_CODE", "density_band"]
    wide = use.pivot_table(
        index=id_cols,
        columns="year",
        values="radiance_analysis",
        aggfunc="first",
        observed=True,
    ).reset_index()

    for year in YEARS:
        if year not in wide.columns:
            wide[year] = np.nan

    wide = wide.rename(
        columns={
            2020: "radiance_2020",
            2021: "radiance_2021",
            2022: "radiance_2022",
        }
    )
    wide["flank_mean"] = 0.5 * (wide["radiance_2020"] + wide["radiance_2022"])
    wide["delta_nw"] = wide["radiance_2021"] - wide["flank_mean"]

    positive = (wide["radiance_2021"] > 0) & (wide["flank_mean"] > 0)
    wide["ratio_2021_to_flanks"] = np.where(
        positive, wide["radiance_2021"] / wide["flank_mean"], np.nan
    )
    wide["pct_change"] = 100.0 * (wide["ratio_2021_to_flanks"] - 1.0)
    wide["increasing"] = wide["delta_nw"] > 0
    wide["paired_valid"] = (
        wide[["radiance_2020", "radiance_2021", "radiance_2022"]]
        .notna()
        .all(axis=1)
    )
    wide["density_band"] = pd.Categorical(
        wide["density_band"], categories=list(DENSITY_ORDER), ordered=True
    )
    return wide


def _summarise_group(group: pd.DataFrame) -> pd.Series:
    paired = group[group["paired_valid"]].copy()
    ratios = paired["ratio_2021_to_flanks"]

    return pd.Series(
        {
            "n_nuts3": paired["NUTS_ID"].nunique(),
            "n_countries": paired["CNTR_CODE"].nunique(),
            "mean_delta_nw": paired["delta_nw"].mean(),
            "median_delta_nw": paired["delta_nw"].median(),
            "mean_pct_change": paired["pct_change"].mean(),
            "median_pct_change": paired["pct_change"].median(),
            "geometric_pct_change": _geometric_pct(ratios),
            "pct_regions_increasing": (
                100.0 * paired["increasing"].mean() if len(paired) else np.nan
            ),
        }
    )


def make_summary_tables(
    region: pd.DataFrame, processed_out: Path
) -> dict[str, pd.DataFrame]:
    overall = (
        region.groupby(["product", "month"], observed=True)
        .apply(_summarise_group)
        .reset_index()
    )

    country_medians = (
        region[region["paired_valid"]]
        .groupby(["product", "month", "CNTR_CODE"], observed=True)["pct_change"]
        .median()
        .rename("country_median_pct")
        .reset_index()
    )
    cb = (
        country_medians.groupby(["product", "month"], observed=True)[
            "country_median_pct"
        ]
        .median()
        .rename("country_balanced_median_pct_change")
        .reset_index()
    )
    overall = overall.merge(cb, on=["product", "month"], how="left")

    by_band = (
        region.groupby(["product", "month", "density_band"], observed=True)
        .apply(_summarise_group)
        .reset_index()
    )

    band_country_medians = (
        region[region["paired_valid"]]
        .groupby(
            ["product", "month", "density_band", "CNTR_CODE"], observed=True
        )["pct_change"]
        .median()
        .rename("country_median_pct")
        .reset_index()
    )
    band_cb = (
        band_country_medians.groupby(
            ["product", "month", "density_band"], observed=True
        )["country_median_pct"]
        .median()
        .rename("country_balanced_median_pct_change")
        .reset_index()
    )
    by_band = by_band.merge(
        band_cb, on=["product", "month", "density_band"], how="left"
    )

    by_country = (
        region.groupby(["product", "month", "CNTR_CODE"], observed=True)
        .apply(_summarise_group)
        .reset_index()
    )

    for frame in (overall, by_band, by_country):
        frame.insert(
            frame.columns.get_loc("month") + 1,
            "month_name",
            frame["month"].map(month_name),
        )

    contrib_frames: list[pd.DataFrame] = []

    for product in PRODUCTS:
        rp = region[
            (region["product"] == product) & region["paired_valid"]
        ].copy()

        complete_counts = rp.groupby("NUTS_ID")["month"].nunique()
        complete_ids = complete_counts[complete_counts == 12].index
        rc = rp[rp["NUTS_ID"].isin(complete_ids)].copy()

        groups: list[tuple[str, pd.DataFrame]] = [("ALL", rc)]
        for band in DENSITY_ORDER:
            groups.append(
                (band, rc[rc["density_band"].astype(str) == band])
            )

        for label, group in groups:
            if group.empty:
                continue

            month_stats = (
                group.groupby("month", observed=True)
                .agg(
                    n_complete_nuts3=("NUTS_ID", "nunique"),
                    n_countries=("CNTR_CODE", "nunique"),
                    mean_delta_nw=("delta_nw", "mean"),
                    median_delta_nw=("delta_nw", "median"),
                    median_pct_change=("pct_change", "median"),
                )
                .reset_index()
            )

            month_stats["product"] = product
            month_stats["density_band"] = label
            month_stats["monthly_contribution_to_mean_nw"] = (
                month_stats["mean_delta_nw"] / 12.0
            )
            net = month_stats["monthly_contribution_to_mean_nw"].sum()
            month_stats["net_reconstructed_mean_delta_nw"] = net

            if np.isfinite(net) and not math.isclose(
                net, 0.0, abs_tol=1e-15
            ):
                month_stats["share_of_net_reconstructed_delta_pct"] = (
                    100.0
                    * month_stats["monthly_contribution_to_mean_nw"]
                    / net
                )
            else:
                month_stats["share_of_net_reconstructed_delta_pct"] = np.nan

            month_stats["cumulative_reconstructed_delta_nw"] = month_stats[
                "monthly_contribution_to_mean_nw"
            ].cumsum()

            contrib_frames.append(month_stats)

    contributions = (
        pd.concat(contrib_frames, ignore_index=True)
        if contrib_frames
        else pd.DataFrame()
    )

    if not contributions.empty:
        contributions.insert(
            contributions.columns.get_loc("month") + 1,
            "month_name",
            contributions["month"].map(month_name),
        )

    region_path = processed_out / "monthly_region_2021_vs_flanks.csv"
    overall_path = processed_out / "monthly_overall_2021_vs_flanks.csv"
    band_path = processed_out / "monthly_density_band_2021_vs_flanks.csv"
    country_path = processed_out / "monthly_country_2021_vs_flanks.csv"
    contrib_path = processed_out / "monthly_complete_case_contributions.csv"

    region.to_csv(region_path, index=False)
    overall.to_csv(overall_path, index=False)
    by_band.to_csv(band_path, index=False)
    by_country.to_csv(country_path, index=False)
    contributions.to_csv(contrib_path, index=False)

    return {
        "region": region,
        "overall": overall,
        "band": by_band,
        "country": by_country,
        "contribution": contributions,
    }


def _display_grid(
    bounds: tuple[float, float, float, float], max_pixels_across: int
):
    west, south, east, north = bounds
    width = int(max_pixels_across)
    height = max(
        1, int(round(width * (north - south) / (east - west)))
    )
    transform = from_bounds(west, south, east, north, width, height)
    return width, height, transform


def read_display_pair(
    pair: RasterPair,
    bounds: tuple[float, float, float, float],
    max_pixels_across: int,
) -> tuple[np.ndarray, np.ndarray]:
    width, height, transform = _display_grid(bounds, max_pixels_across)
    rad = np.full((height, width), np.nan, dtype=np.float32)
    cov = np.zeros((height, width), dtype=np.float32)

    with rasterio.open(pair.radiance) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=rad,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=transform,
            dst_crs="EPSG:4326",
            dst_nodata=np.nan,
            resampling=Resampling.nearest,
        )

    with rasterio.open(pair.cloudfree) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=cov,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=transform,
            dst_crs="EPSG:4326",
            dst_nodata=0,
            resampling=Resampling.nearest,
        )

    valid = np.isfinite(rad) & np.isfinite(cov) & (cov > 0)
    rad[~valid] = np.nan
    return rad, valid


def _robust_high(
    arrays: Iterable[np.ndarray], q: float, *, absolute: bool = False
) -> float:
    vals: list[float] = []

    for arr in arrays:
        x = arr[np.isfinite(arr)]
        if absolute:
            x = np.abs(x)
        else:
            x = x[x > 0]

        if x.size:
            vals.append(float(np.quantile(x, q)))

    if not vals:
        return 1.0

    high = max(vals)
    return high if high > 0 else 1.0


def _plot_boundaries(ax, countries: gpd.GeoDataFrame) -> None:
    countries.boundary.plot(
        ax=ax,
        color="0.75",
        linewidth=0.25,
        alpha=0.7,
        zorder=3,
    )
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_raw_year_atlas(
    product: str,
    year: int,
    arrays: dict[tuple[int, int], np.ndarray],
    countries: gpd.GeoDataFrame,
    bounds: tuple[float, float, float, float],
    vmax: float,
    output: Path,
) -> None:
    norm = SymLogNorm(
        linthresh=0.2,
        linscale=0.8,
        vmin=0.0,
        vmax=vmax,
        base=10,
    )

    fig, axes = plt.subplots(
        3, 4, figsize=(15.2, 10.8), constrained_layout=True
    )
    axes = axes.ravel()
    image = None

    for month, ax in enumerate(axes, start=1):
        arr = arrays[(year, month)].copy()
        arr[arr < 0] = 0.0

        image = ax.imshow(
            arr,
            extent=(bounds[0], bounds[2], bounds[1], bounds[3]),
            origin="upper",
            cmap="magma",
            norm=norm,
            interpolation="nearest",
            rasterized=True,
        )
        _plot_boundaries(ax, countries)
        ax.set_title(month_name(month), fontsize=10)

    cbar = fig.colorbar(
        image,
        ax=axes.tolist(),
        orientation="horizontal",
        shrink=0.70,
        pad=0.02,
    )
    cbar.set_label(
        r"Monthly average radiance (nW cm$^{-2}$ sr$^{-1}$)"
    )

    fig.suptitle(
        f"{PRODUCT_LABELS[product]} — {year}",
        fontsize=15,
        weight="semibold",
    )
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_monthly_anomaly_atlas(
    product: str,
    arrays: dict[tuple[int, int], np.ndarray],
    countries: gpd.GeoDataFrame,
    bounds: tuple[float, float, float, float],
    output_delta: Path,
    output_pct: Path,
    pct_min_radiance: float,
) -> None:
    deltas: dict[int, np.ndarray] = {}
    pcts: dict[int, np.ndarray] = {}

    for month in range(1, 13):
        a20 = arrays[(2020, month)]
        a21 = arrays[(2021, month)]
        a22 = arrays[(2022, month)]

        valid = np.isfinite(a20) & np.isfinite(a21) & np.isfinite(a22)

        flank = np.full_like(a21, np.nan)
        flank[valid] = 0.5 * (a20[valid] + a22[valid])

        delta = np.full_like(a21, np.nan)
        delta[valid] = a21[valid] - flank[valid]
        deltas[month] = delta

        pct = np.full_like(a21, np.nan)
        pct_valid = (
            valid
            & (flank >= pct_min_radiance)
            & (a21 >= 0)
        )
        pct[pct_valid] = (
            100.0 * delta[pct_valid] / flank[pct_valid]
        )
        pcts[month] = pct

    delta_limit = _robust_high(
        deltas.values(), 0.99, absolute=True
    )
    pct_limit = _robust_high(
        pcts.values(), 0.985, absolute=True
    )
    pct_limit = min(max(pct_limit, 10.0), 200.0)

    plot_specs = [
        (
            deltas,
            delta_limit,
            output_delta,
            (
                r"2021 − mean(2020, 2022) radiance "
                r"(nW cm$^{-2}$ sr$^{-1}$)"
            ),
            "absolute anomaly",
        ),
        (
            pcts,
            pct_limit,
            output_pct,
            (
                "2021 vs flanks (%) — displayed where "
                f"flank radiance ≥ {pct_min_radiance:g}"
            ),
            "percentage anomaly",
        ),
    ]

    for data, limit, output, label, title_suffix in plot_specs:
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)

        fig, axes = plt.subplots(
            3, 4, figsize=(15.2, 10.8), constrained_layout=True
        )
        axes = axes.ravel()
        image = None

        for month, ax in enumerate(axes, start=1):
            image = ax.imshow(
                data[month],
                extent=(bounds[0], bounds[2], bounds[1], bounds[3]),
                origin="upper",
                cmap="coolwarm",
                norm=norm,
                interpolation="nearest",
                rasterized=True,
            )
            _plot_boundaries(ax, countries)
            ax.set_title(month_name(month), fontsize=10)

        cbar = fig.colorbar(
            image,
            ax=axes.tolist(),
            orientation="horizontal",
            shrink=0.70,
            pad=0.02,
        )
        cbar.set_label(label)

        fig.suptitle(
            (
                f"{PRODUCT_LABELS[product]} — 2021 {title_suffix} "
                f"vs same-month 2020/2022 flanks"
            ),
            fontsize=14,
            weight="semibold",
        )
        fig.savefig(output, dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_product_comparison_2021(
    vcm: dict[tuple[int, int], np.ndarray],
    vcmsl: dict[tuple[int, int], np.ndarray],
    vcm_valid: dict[tuple[int, int], np.ndarray],
    vcmsl_valid: dict[tuple[int, int], np.ndarray],
    countries: gpd.GeoDataFrame,
    bounds: tuple[float, float, float, float],
    diff_output: Path,
    coverage_output: Path,
) -> None:
    diffs: dict[int, np.ndarray] = {}
    coverdiff: dict[int, np.ndarray] = {}

    for month in range(1, 13):
        a = vcm[(2021, month)]
        b = vcmsl[(2021, month)]

        common = np.isfinite(a) & np.isfinite(b)
        d = np.full_like(a, np.nan)
        d[common] = b[common] - a[common]
        diffs[month] = d

        cv = (
            vcmsl_valid[(2021, month)].astype(np.int8)
            - vcm_valid[(2021, month)].astype(np.int8)
        )
        coverdiff[month] = cv.astype(float)

    limit = _robust_high(diffs.values(), 0.99, absolute=True)
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)

    fig, axes = plt.subplots(
        3, 4, figsize=(15.2, 10.8), constrained_layout=True
    )
    axes = axes.ravel()
    image = None

    for month, ax in enumerate(axes, start=1):
        image = ax.imshow(
            diffs[month],
            extent=(bounds[0], bounds[2], bounds[1], bounds[3]),
            origin="upper",
            cmap="coolwarm",
            norm=norm,
            interpolation="nearest",
            rasterized=True,
        )
        _plot_boundaries(ax, countries)
        ax.set_title(month_name(month), fontsize=10)

    cbar = fig.colorbar(
        image,
        ax=axes.tolist(),
        orientation="horizontal",
        shrink=0.70,
        pad=0.02,
    )
    cbar.set_label(
        r"VCMSLCFG − VCMCFG radiance on common-valid pixels "
        r"(nW cm$^{-2}$ sr$^{-1}$)"
    )
    fig.suptitle(
        "2021 monthly product difference on common valid coverage",
        fontsize=14,
        weight="semibold",
    )
    fig.savefig(diff_output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(
        3, 4, figsize=(15.2, 10.8), constrained_layout=True
    )
    axes = axes.ravel()
    image = None

    for month, ax in enumerate(axes, start=1):
        image = ax.imshow(
            coverdiff[month],
            extent=(bounds[0], bounds[2], bounds[1], bounds[3]),
            origin="upper",
            cmap="coolwarm",
            vmin=-1,
            vmax=1,
            interpolation="nearest",
            rasterized=True,
        )
        _plot_boundaries(ax, countries)
        ax.set_title(month_name(month), fontsize=10)

    cbar = fig.colorbar(
        image,
        ax=axes.tolist(),
        orientation="horizontal",
        shrink=0.70,
        pad=0.02,
        ticks=[-1, 0, 1],
    )
    cbar.ax.set_xticklabels(
        ["VCMCFG only", "same", "VCMSLCFG only"]
    )
    fig.suptitle(
        "2021 monthly valid-coverage difference between products",
        fontsize=14,
        weight="semibold",
    )
    fig.savefig(coverage_output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_maps(
    pairs: dict[tuple[str, int, int], RasterPair],
    nuts: gpd.GeoDataFrame,
    figure_out: Path,
    bounds: tuple[float, float, float, float],
    max_pixels_across: int,
    pct_min_radiance: float,
) -> None:
    figure_out.mkdir(parents=True, exist_ok=True)

    west, south, east, north = bounds
    map_box = gpd.GeoDataFrame(
        geometry=[box(west, south, east, north)],
        crs="EPSG:4326",
    )
    europe_nuts = gpd.clip(nuts, map_box)
    countries = europe_nuts.dissolve(by="CNTR_CODE")

    product_arrays: dict[
        str, dict[tuple[int, int], np.ndarray]
    ] = {}
    product_valid: dict[
        str, dict[tuple[int, int], np.ndarray]
    ] = {}

    for product in PRODUCTS:
        arrays: dict[tuple[int, int], np.ndarray] = {}
        valids: dict[tuple[int, int], np.ndarray] = {}

        print(f"Preparing display rasters for {product}...")

        for year in YEARS:
            for month in range(1, 13):
                arr, valid = read_display_pair(
                    pairs[(product, year, month)],
                    bounds,
                    max_pixels_across,
                )
                arrays[(year, month)] = arr
                valids[(year, month)] = valid

        product_arrays[product] = arrays
        product_valid[product] = valids

        vmax = _robust_high(
            arrays.values(), 0.997, absolute=False
        )

        for year in YEARS:
            output = figure_out / f"{product}_raw_monthly_{year}.pdf"
            plot_raw_year_atlas(
                product,
                year,
                arrays,
                countries,
                bounds,
                vmax,
                output,
            )
            print(f"Wrote: {output}")

        delta_output = (
            figure_out
            / f"{product}_2021_minus_2020_2022_flanks.pdf"
        )
        pct_output = (
            figure_out
            / f"{product}_2021_pct_vs_2020_2022_flanks.pdf"
        )

        plot_monthly_anomaly_atlas(
            product,
            arrays,
            countries,
            bounds,
            delta_output,
            pct_output,
            pct_min_radiance,
        )
        print(f"Wrote: {delta_output}")
        print(f"Wrote: {pct_output}")

    diff_output = (
        figure_out
        / "2021_vcmslcfg_minus_vcmcfg_common_coverage.pdf"
    )
    coverage_output = (
        figure_out
        / "2021_vcmslcfg_vs_vcmcfg_coverage_difference.pdf"
    )

    plot_product_comparison_2021(
        product_arrays["vcmcfg"],
        product_arrays["vcmslcfg"],
        product_valid["vcmcfg"],
        product_valid["vcmslcfg"],
        countries,
        bounds,
        diff_output,
        coverage_output,
    )
    print(f"Wrote: {diff_output}")
    print(f"Wrote: {coverage_output}")


def write_summary_text(
    tables: dict[str, pd.DataFrame],
    monthly: pd.DataFrame,
    output: Path,
    min_valid_fraction: float,
) -> None:
    overall = tables["overall"].copy()
    band = tables["band"].copy()
    contribution = tables["contribution"].copy()

    lines: list[str] = []
    lines.append(
        "MONTHLY VIIRS DIAGNOSTIC OF THE 2021 ANNUAL ANOMALY"
    )
    lines.append("=" * 78)
    lines.append("")
    lines.append("LOCKED DESIGN")
    lines.append("- fixed geography: NUTS2016")
    lines.append(
        f"- locked cohort: {EXPECTED_COHORT_SIZE:,} NUTS3 regions"
    )
    lines.append(
        "- comparison: 2021 month vs mean of the same month in 2020 and 2022"
    )
    lines.append("- products: VCMCFG and VCMSLCFG")
    lines.append(
        "- cf_cvg > 0 masks pixels without cloud-free observations"
    )
    lines.append(
        "- inferential region-month summaries require >= "
        f"{100 * min_valid_fraction:.0f}% valid NUTS3 area"
    )
    lines.append(
        "- zonal means use exact polygon overlap with spherical physical-area weighting"
    )
    lines.append("")

    lines.append("VALID REGION-MONTH SUPPORT")
    support = (
        monthly.groupby(["product", "year", "month"], observed=True)[
            "analysis_valid"
        ]
        .sum()
        .astype(int)
        .reset_index(name="n_valid_nuts3")
    )
    lines.append(support.to_string(index=False))
    lines.append("")

    lines.append("OVERALL MONTHLY 2021 VS FLANKS")
    lines.append(
        overall.to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    lines.append("")

    lines.append("BY FIXED POPULATION-DENSITY BAND")
    lines.append(
        band.to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    lines.append("")

    if not contribution.empty:
        lines.append("COMPLETE-CASE MONTHLY CONTRIBUTION DIAGNOSTIC")
        lines.append(
            "The additive contribution table reconstructs the change in the mean "
            "of 12 monthly regional means among regions with all 12 paired months. "
            "It is NOT the EOG annual compositing algorithm and must not be used "
            "as a replacement annual product."
        )
        lines.append(
            contribution.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        lines.append("")

    lines.append("INTERPRETATION RULE")
    lines.append(
        "If the 2021 anomaly is concentrated in one or a few winter months and "
        "in coherent geographic areas, especially where VCMSLCFG adds "
        "high-latitude coverage, treat the annual pulse as a surface-state/"
        "coverage/product diagnostic before fitting a biological or long-term "
        "temporal break. If it is broad across months and both products, the "
        "annual pulse is less likely to be explained by a small seasonal subset."
    )
    lines.append("")
    lines.append(
        "This stage is diagnostic. It does not overwrite annual VIIRS, remove "
        "2021, recalibrate the production time series, or alter the locked "
        "NUTS2016 cohort."
    )

    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="VIIRS raw-data root scanned recursively. Default: <viirs_raw>.",
    )
    parser.add_argument(
        "--min-valid-area-fraction",
        type=float,
        default=DEFAULT_MIN_VALID_AREA_FRACTION,
        help=(
            "Minimum NUTS3 area with cf_cvg > 0 required for "
            "region-month summaries."
        ),
    )
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--skip-zonal", action="store_true")
    parser.add_argument("--skip-maps", action="store_true")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute zonal checkpoints.",
    )
    parser.add_argument(
        "--max-pixels-across",
        type=int,
        default=DEFAULT_MAX_PIXELS_ACROSS,
    )
    parser.add_argument(
        "--percent-map-min-radiance",
        type=float,
        default=DEFAULT_PERCENT_MAP_MIN_RADIANCE,
        help=(
            "Display-only minimum flank radiance for percentage anomaly maps."
        ),
    )
    args = parser.parse_args()

    if not (0.0 < args.min_valid_area_fraction <= 1.0):
        raise ValueError(
            "--min-valid-area-fraction must be in (0, 1]."
        )

    config = load_config(args.config)
    viirs_raw = configured_path(config, "viirs_raw")
    processed = configured_path(config, "processed")
    figures = configured_path(config, "figures")

    input_dir = args.input_dir or viirs_raw
    if not input_dir.is_absolute():
        input_dir = Path.cwd() / input_dir

    processed_out = processed / "fixed_nuts2016_monthly_product_atlas"
    figure_out = (
        figures / "population_relationship" / "fixed_nuts2016_monthly_product_atlas"
    )
    processed_out.mkdir(parents=True, exist_ok=True)
    figure_out.mkdir(parents=True, exist_ok=True)

    pairs, missing = discover_inputs(input_dir)
    print_inventory(pairs, missing, input_dir)
    if args.inventory_only:
        return

    if missing:
        raise FileNotFoundError(
            "Monthly input inventory is incomplete. Download/extract the missing "
            "EOG monthly average-radiance and cf_cvg GeoTIFFs, then rerun."
        )

    # Reuse the canonical locked fixed-NUTS2016 panel already established by
    # the future-radiance workflow.  Do not reconstruct or silently change the
    # 1,514-region cohort or its fixed density classes here.
    panel = load_locked_panel(config)
    cohort = (
        panel.loc[panel["year"].astype(int) == 2020, ["NUTS_ID", "density_band"]]
        .drop_duplicates("NUTS_ID")
        .copy()
    )
    if len(cohort) != EXPECTED_COHORT_SIZE:
        raise RuntimeError(
            f"Expected {EXPECTED_COHORT_SIZE:,} locked NUTS3 regions, "
            f"found {len(cohort):,}."
        )

    release = int(
        config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"]
    )
    nuts = load_fixed_nuts(config, release).copy()
    id_field = config["nuts"]["id_field"]
    if id_field != "NUTS_ID":
        nuts = nuts.rename(columns={id_field: "NUTS_ID"})
    nuts["NUTS_ID"] = nuts["NUTS_ID"].astype(str)
    keep = [
        x for x in ["NUTS_ID", "CNTR_CODE", "NUTS_NAME", "NAME_LATN", "geometry"]
        if x in nuts.columns
    ]
    nuts = nuts[keep].merge(
        cohort, on="NUTS_ID", how="inner", validate="one_to_one"
    )
    if len(nuts) != EXPECTED_COHORT_SIZE:
        raise RuntimeError(
            f"NUTS2016 geometry contains {len(nuts):,}/"
            f"{EXPECTED_COHORT_SIZE:,} locked regions."
        )
    if "CNTR_CODE" not in nuts.columns:
        nuts["CNTR_CODE"] = nuts["NUTS_ID"].str[:2]

    map_bounds = tuple(
        float(v)
        for v in config.get("analysis", {})
        .get("europe_maps", {})
        .get("extent_4326", DEFAULT_MAP_BOUNDS)
    )

    combined_path = (
        processed_out / "monthly_nuts2016_radiance_2020_2022.parquet"
    )

    if not args.skip_zonal:
        monthly = run_zonal_extraction(
            pairs,
            nuts,
            processed_out,
            args.min_valid_area_fraction,
            args.force,
        )
        region = build_region_anomalies(monthly)
        tables = make_summary_tables(region, processed_out)
        write_summary_text(
            tables,
            monthly,
            processed_out / "monthly_2021_diagnostic_summary.txt",
            args.min_valid_area_fraction,
        )

    elif combined_path.exists():
        monthly = pd.read_parquet(combined_path)
        region = build_region_anomalies(monthly)
        tables = make_summary_tables(region, processed_out)
        write_summary_text(
            tables,
            monthly,
            processed_out / "monthly_2021_diagnostic_summary.txt",
            args.min_valid_area_fraction,
        )

    if not args.skip_maps:
        run_maps(
            pairs,
            nuts,
            figure_out,
            map_bounds,
            args.max_pixels_across,
            args.percent_map_min_radiance,
        )

    print("\nDONE")
    print(f"Processed outputs: {processed_out}")
    print(f"Figures: {figure_out}")


if __name__ == "__main__":
    main()
