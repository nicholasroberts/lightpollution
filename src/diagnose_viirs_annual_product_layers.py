#!/usr/bin/env python3
"""Audit VIIRS annual product layers around the 2021 pulse.

This is a pre-specified product-quality diagnostic, not a new biological trend
model. It uses the same locked NUTS2016 cohort and compares three annual VNL
layers for 2020-2023:

- average_masked: the existing raw radiance response already extracted;
- median_masked: a robust annual radiance summary;
- cf_cvg: number of cloud-free observations contributing to each pixel.

The diagnostic tests whether the 2021 average-masked pulse is also present in
median-masked radiance and whether region-level radiance anomalies are related
to changes in cloud-free observation coverage. No production file is changed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
from scipy.stats import linregress, spearmanr

from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import (
    load_fixed_nuts,
    prepare_analysis_panel,
)


YEARS = [2020, 2021, 2022, 2023]
LAYERS = ["median_masked", "cf_cvg"]


def restore_density_order(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    panel = panel.copy()
    panel["density_band"] = pd.Categorical(
        panel["density_band"], categories=labels, ordered=True
    )
    if panel["density_band"].isna().any():
        raise RuntimeError("Unexpected density-band label in locked panel.")
    return panel


def load_locked_panel(config: dict) -> pd.DataFrame:
    cfg = config["analysis"]["fixed_nuts2016_2014_2024"]
    diag = config["analysis"]["fixed_nuts2016_2021_diagnostic"]
    processed = configured_path(config, "processed")
    fixed = pd.read_parquet(processed / diag["fixed_radiance_file"]).copy()
    reference = pd.read_csv(processed / diag["population_reference_file"]).copy()
    panel = prepare_analysis_panel(fixed, reference, config)
    panel = restore_density_order(panel, config)
    panel = panel[panel["year"].isin(YEARS)].copy()

    expected_n = int(diag["expected_n_nuts3"])
    counts = panel.groupby("NUTS_ID")["year"].nunique()
    keep = counts[counts == len(YEARS)].index
    panel = panel[panel["NUTS_ID"].isin(keep)].copy()
    if panel["NUTS_ID"].nunique() != expected_n:
        raise RuntimeError(
            f"Product-layer audit expected {expected_n} fixed regions; "
            f"found {panel['NUTS_ID'].nunique()}."
        )
    return panel.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)


def add_latitude(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release)[["NUTS_ID", "geometry"]].copy()
    projected = nuts.to_crs("EPSG:3035")
    centres = gpd.GeoDataFrame(
        {"NUTS_ID": projected["NUTS_ID"].values},
        geometry=projected.geometry.centroid,
        crs=projected.crs,
    ).to_crs("EPSG:4326")
    loc = pd.DataFrame(
        {
            "NUTS_ID": centres["NUTS_ID"].values,
            "centroid_latitude": centres.geometry.y.to_numpy(float),
        }
    )
    return panel.merge(loc, on="NUTS_ID", how="left", validate="many_to_one")


def expected_companion_name(avg_name: str, layer: str) -> str:
    suffix = ".average_masked.dat.tif"
    if not avg_name.endswith(suffix):
        raise RuntimeError(
            f"Cannot derive companion layer from unexpected source name: {avg_name}"
        )
    return avg_name[: -len(suffix)] + f".{layer}.dat.tif"


def _fallback_companion_matches(
    root: Path,
    year: int,
    layer: str,
) -> list[Path]:
    """Find a decompressed companion despite harmless filename variations.

    Browser downloads and EOG re-packaging can preserve the scientific layer
    while changing the creation timestamp or adding a duplicate-download
    suffix. The observation year and layer token are therefore the scientific
    identifiers used for the fallback search. Exact expected filenames remain
    preferred.
    """
    year_token = f"_{year}_global_"
    layer_tokens = {
        f".{layer}.",
        f".{layer.replace('_', '-')}.",
    }
    matches = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        name = path.name
        if not name.lower().endswith((".tif", ".tiff")):
            continue
        if year_token not in name:
            continue
        if not any(token in name for token in layer_tokens):
            continue
        matches.append(path)
    return sorted(matches)


def locate_companions(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    root = configured_path(config, "viirs_raw")
    rows = []
    missing = []
    for year in YEARS:
        g = panel[panel["year"] == year]
        avg_names = sorted(set(g["viirs_file"].astype(str)))
        if len(avg_names) != 1:
            raise RuntimeError(
                f"Expected one average_masked source for {year}; found {avg_names}"
            )
        avg_name = avg_names[0]
        for layer in LAYERS:
            expected = expected_companion_name(avg_name, layer)
            exact_matches = list(root.rglob(expected))
            match_mode = "exact"

            if len(exact_matches) > 1:
                raise RuntimeError(
                    f"Multiple exact files match {expected}: "
                    + ", ".join(str(p) for p in exact_matches)
                )

            if exact_matches:
                matches = exact_matches
            else:
                matches = _fallback_companion_matches(root, year, layer)
                match_mode = "year_layer_fallback"

            if len(matches) > 1:
                raise RuntimeError(
                    f"Multiple decompressed candidate files found for "
                    f"{year} {layer}. Keep or rename the intended file so the "
                    "choice is unambiguous:\n  "
                    + "\n  ".join(str(p) for p in matches)
                )

            if not matches:
                missing.append(expected)
                path = None
                match_mode = "missing"
            else:
                path = matches[0]

            rows.append(
                {
                    "year": year,
                    "layer": layer,
                    "expected_filename": expected,
                    "actual_filename": path.name if path is not None else "",
                    "match_mode": match_mode,
                    "found": path is not None,
                    "relative_path": (
                        str(path.relative_to(root)) if path is not None else ""
                    ),
                    "size_bytes": (
                        int(path.stat().st_size) if path is not None else np.nan
                    ),
                }
            )

    manifest = pd.DataFrame(rows)
    if missing:
        needed = "\n".join(f"  {name}" for name in missing)
        raise FileNotFoundError(
            "Missing required annual VNL companion rasters under "
            f"{root}:\n{needed}\n"
            "No unique decompressed file with the same observation year and "
            "layer token was found either. Check the local filenames with:\n"
            "  find data/raw/viirs -maxdepth 2 -type f "
            "\\( -iname '*median*' -o -iname '*cf*cvg*' \\) -printf '%f\\n'\n"
        )
    return manifest


def companion_paths(manifest: pd.DataFrame, config: dict) -> dict[tuple[int, str], Path]:
    root = configured_path(config, "viirs_raw")
    out = {}
    for _, row in manifest.iterrows():
        out[(int(row["year"]), str(row["layer"]))] = (
            root / str(row["relative_path"])
        )
    return out


def extract_layer(
    path: Path,
    layer: str,
    year: int,
    nuts: gpd.GeoDataFrame,
    keep_ids: set[str],
    config: dict,
) -> pd.DataFrame:
    print(f"Extracting {year} {layer}: {path.name}", flush=True)
    id_field = config["nuts"]["id_field"]
    with rasterio.open(path) as ds:
        if ds.crs is None or ds.crs.to_string() != config["nuts"]["crs"]:
            raise RuntimeError(f"Unexpected CRS in {path}")
        source = RasterioRasterSource(ds, 1, name=layer)
        result = exact_extract(
            source,
            nuts,
            ["mean(coverage_weight=area_spherical_m2)"],
            include_cols=[id_field],
            output="pandas",
            progress=True,
        )
    mean_col = f"{layer}_mean"
    # exactextract prefixes statistics with the RasterSource name only when
    # multiple raster sources are supplied in the same call. This diagnostic
    # intentionally extracts one companion layer at a time, for which the
    # returned statistic is simply named "mean". Normalize that single-source
    # output to an explicit layer-specific column here.
    if mean_col not in result.columns:
        if "mean" in result.columns:
            result = result.rename(columns={"mean": mean_col})
        else:
            raise RuntimeError(
                f"Expected 'mean' or {mean_col} from exactextract; "
                f"found {list(result.columns)}"
            )
    result[id_field] = result[id_field].astype(str)
    result = result[result[id_field].isin(keep_ids)][[id_field, mean_col]].copy()
    if result[id_field].duplicated().any():
        raise RuntimeError(f"Duplicate {id_field} rows in {year} {layer} extraction.")
    observed_ids = set(result[id_field])
    if observed_ids != keep_ids:
        missing_ids = sorted(keep_ids - observed_ids)
        extra_ids = sorted(observed_ids - keep_ids)
        raise RuntimeError(
            f"{year} {layer} extraction does not match the locked cohort: "
            f"missing={len(missing_ids)}, extra={len(extra_ids)}."
        )
    result.insert(0, "year", year)
    return result


def _validate_cached_layer(
    cached: pd.DataFrame,
    year: int,
    layer: str,
    keep_ids: set[str],
) -> pd.DataFrame:
    mean_col = f"{layer}_mean"
    required = {"year", "NUTS_ID", mean_col}
    missing_cols = required - set(cached.columns)
    if missing_cols:
        raise RuntimeError(
            f"Cached {year} {layer} is missing columns {sorted(missing_cols)}."
        )
    out = cached[["year", "NUTS_ID", mean_col]].copy()
    out["NUTS_ID"] = out["NUTS_ID"].astype(str)
    if set(out["year"].astype(int)) != {int(year)}:
        raise RuntimeError(f"Cached {layer} file contains the wrong year.")
    if out["NUTS_ID"].duplicated().any():
        raise RuntimeError(f"Cached {year} {layer} has duplicate NUTS_ID rows.")
    if set(out["NUTS_ID"]) != keep_ids:
        raise RuntimeError(
            f"Cached {year} {layer} does not match the locked NUTS3 cohort."
        )
    return out


def extract_all(
    panel: pd.DataFrame,
    manifest: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release)
    keep_ids = set(panel["NUTS_ID"].astype(str))
    paths = companion_paths(manifest, config)

    # Cache every expensive year-layer extraction immediately. Downstream
    # diagnostics are cheap; a later coding error must never force the user to
    # repeat all eight global-raster extractions.
    cache_dir = configured_path(config, "processed") / "viirs_product_layer_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    pieces = []
    for year in YEARS:
        year_parts = []
        for layer in LAYERS:
            cache_path = cache_dir / f"nuts2016_{year}_{layer}.parquet"
            if cache_path.exists():
                print(f"Using cached {year} {layer}: {cache_path}", flush=True)
                part = _validate_cached_layer(
                    pd.read_parquet(cache_path), year, layer, keep_ids
                )
            else:
                part = extract_layer(
                    paths[(year, layer)],
                    layer,
                    year,
                    nuts,
                    keep_ids,
                    config,
                )
                part.to_parquet(cache_path, index=False)
                print(f"Cached {year} {layer}: {cache_path}", flush=True)
            year_parts.append(part)

        merged = year_parts[0]
        for part in year_parts[1:]:
            merged = merged.merge(
                part, on=["year", "NUTS_ID"], how="inner", validate="one_to_one"
            )
        if len(merged) != len(keep_ids):
            raise RuntimeError(
                f"Merged companion extraction for {year} has {len(merged)} "
                f"rows; expected {len(keep_ids)}."
            )
        pieces.append(merged)

    out = pd.concat(pieces, ignore_index=True)
    expected_rows = len(keep_ids) * len(YEARS)
    if len(out) != expected_rows:
        raise RuntimeError(
            f"Companion table has {len(out)} rows; expected {expected_rows}."
        )
    return out


def build_analysis(panel: pd.DataFrame, companions: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "NUTS_ID",
        "CNTR_CODE",
        "density_band",
        "reference_population_density",
        "centroid_latitude",
        "year",
        "radiance_mean_raw",
    ]
    d = panel[cols].merge(
        companions,
        on=["year", "NUTS_ID"],
        how="left",
        validate="one_to_one",
    )
    for col in ["radiance_mean_raw", "median_masked_mean", "cf_cvg_mean"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
        if (~np.isfinite(d[col])).any():
            raise RuntimeError(f"Non-finite values in {col}.")

    # Radiance layers must be positive for log-change analysis; the locked
    # average-masked cohort already enforces that. For median-masked, record
    # zero values rather than silently dropping them.
    d["median_positive"] = d["median_masked_mean"] > 0
    return d.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)


def make_transition_table(d: pd.DataFrame) -> pd.DataFrame:
    x = d.copy()
    for metric, col in [
        ("average_masked", "radiance_mean_raw"),
        ("median_masked", "median_masked_mean"),
    ]:
        x[f"prev_{metric}"] = x.groupby("NUTS_ID")[col].shift(1)
        valid = (x[col] > 0) & (x[f"prev_{metric}"] > 0)
        x[f"delta_log10_{metric}"] = np.nan
        x.loc[valid, f"delta_log10_{metric}"] = np.log10(
            x.loc[valid, col] / x.loc[valid, f"prev_{metric}"]
        )
        x[f"pct_change_{metric}"] = 100.0 * (
            10.0 ** x[f"delta_log10_{metric}"] - 1.0
        )

    x["prev_cf_cvg"] = x.groupby("NUTS_ID")["cf_cvg_mean"].shift(1)
    x["delta_cf_cvg"] = x["cf_cvg_mean"] - x["prev_cf_cvg"]
    x["pct_change_cf_cvg"] = np.nan
    valid_cf = x["prev_cf_cvg"] > 0
    x.loc[valid_cf, "pct_change_cf_cvg"] = 100.0 * (
        x.loc[valid_cf, "cf_cvg_mean"] / x.loc[valid_cf, "prev_cf_cvg"] - 1.0
    )
    x["previous_year"] = x.groupby("NUTS_ID")["year"].shift(1)
    return x[x["previous_year"].notna()].copy()


def geom_change(g: pd.DataFrame, col: str) -> float:
    vals = g[col].dropna()
    return float(100.0 * (10.0 ** vals.mean() - 1.0))


def annual_summary(trans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for year, g in trans.groupby("year", observed=True):
        rows.append(
            {
                "end_year": int(year),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "average_masked_geometric_pct": geom_change(
                    g, "delta_log10_average_masked"
                ),
                "median_masked_geometric_pct": geom_change(
                    g, "delta_log10_median_masked"
                ),
                "median_masked_valid_n": int(
                    g["delta_log10_median_masked"].notna().sum()
                ),
                "mean_cf_cvg": float(g["cf_cvg_mean"].mean()),
                "mean_delta_cf_cvg": float(g["delta_cf_cvg"].mean()),
                "median_delta_cf_cvg": float(g["delta_cf_cvg"].median()),
            }
        )
    return pd.DataFrame(rows).sort_values("end_year")


def density_summary(trans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, band), g in trans.groupby(
        ["year", "density_band"], observed=True
    ):
        rows.append(
            {
                "end_year": int(year),
                "density_band": str(band),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "average_masked_geometric_pct": geom_change(
                    g, "delta_log10_average_masked"
                ),
                "median_masked_geometric_pct": geom_change(
                    g, "delta_log10_median_masked"
                ),
                "mean_delta_cf_cvg": float(g["delta_cf_cvg"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["end_year", "density_band"])


def _spearman(a: pd.Series, b: pd.Series) -> tuple[float, float, int]:
    ok = np.isfinite(a.to_numpy(float)) & np.isfinite(b.to_numpy(float))
    r = spearmanr(a[ok], b[ok])
    return float(r.statistic), float(r.pvalue), int(ok.sum())


def diagnostic_table(trans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for year in [2021, 2022, 2023]:
        g = trans[trans["year"] == year].copy()
        rho_am, p_am, n_am = _spearman(
            g["delta_log10_average_masked"],
            g["delta_log10_median_masked"],
        )
        rho_cf, p_cf, n_cf = _spearman(
            g["delta_log10_average_masked"], g["delta_cf_cvg"]
        )
        rho_lat, p_lat, n_lat = _spearman(
            g["delta_log10_average_masked"], g["centroid_latitude"]
        )
        rows.append(
            {
                "end_year": year,
                "rho_average_vs_median_change": rho_am,
                "p_average_vs_median_change": p_am,
                "n_average_vs_median": n_am,
                "rho_average_change_vs_delta_cf_cvg": rho_cf,
                "p_average_change_vs_delta_cf_cvg": p_cf,
                "n_average_vs_cf_cvg": n_cf,
                "rho_average_change_vs_latitude": rho_lat,
                "p_average_change_vs_latitude": p_lat,
                "n_average_vs_latitude": n_lat,
            }
        )
    return pd.DataFrame(rows)


def net_2020_2022(d: pd.DataFrame) -> pd.DataFrame:
    """Two-year endpoint change without mixing single- and multi-level columns."""
    meta = (
        d[["NUTS_ID", "density_band", "CNTR_CODE"]]
        .drop_duplicates("NUTS_ID")
        .copy()
    )
    if meta["NUTS_ID"].duplicated().any():
        raise RuntimeError("Duplicate NUTS_ID metadata rows.")

    subset = d[d["year"].isin([2020, 2022])][
        ["NUTS_ID", "year", "radiance_mean_raw", "median_masked_mean"]
    ].copy()
    if subset.duplicated(["NUTS_ID", "year"]).any():
        raise RuntimeError("Duplicate NUTS_ID-year rows in 2020/2022 endpoints.")

    wide = subset.pivot(
        index="NUTS_ID",
        columns="year",
        values=["radiance_mean_raw", "median_masked_mean"],
    )
    expected_pairs = [
        ("radiance_mean_raw", 2020),
        ("radiance_mean_raw", 2022),
        ("median_masked_mean", 2020),
        ("median_masked_mean", 2022),
    ]
    missing_pairs = [pair for pair in expected_pairs if pair not in wide.columns]
    if missing_pairs:
        raise RuntimeError(
            f"Missing 2020/2022 endpoint columns: {missing_pairs}"
        )

    # Flatten the pivoted MultiIndex before combining it with ordinary
    # metadata columns. Pandas 3 correctly rejects a direct one-level/two-level
    # join, which was the source of the previous failure.
    wide.columns = [
        f"{metric}_{int(year)}" for metric, year in wide.columns
    ]
    wide = wide.reset_index().merge(
        meta,
        on="NUTS_ID",
        how="left",
        validate="one_to_one",
    )

    rows = []
    for band, g in wide.groupby("density_band", observed=True):
        row = {"density_band": str(band), "n_nuts3": int(len(g))}
        for label, metric in [
            ("average_masked", "radiance_mean_raw"),
            ("median_masked", "median_masked_mean"),
        ]:
            a = g[f"{metric}_2020"]
            b = g[f"{metric}_2022"]
            ok = (a > 0) & (b > 0) & np.isfinite(a) & np.isfinite(b)
            if not ok.any():
                raise RuntimeError(
                    f"No valid 2020/2022 {label} endpoints for density band {band}."
                )
            delta = np.log10(b[ok] / a[ok])
            row[f"net_2020_2022_geometric_pct_{label}"] = float(
                100.0 * (10.0 ** delta.mean() - 1.0)
            )
            row[f"n_valid_{label}"] = int(ok.sum())
        rows.append(row)
    return pd.DataFrame(rows)


def downstream_self_test() -> None:
    """Cheap synthetic smoke test for every post-extraction transformation.

    This runs before any global raster extraction. It is designed to catch
    pandas reshaping/join mistakes and schema drift before expensive work.
    """
    rows = []
    specs = [
        ("A", "AA", "<30", 10.0, 45.0),
        ("B", "AA", "30-100", 60.0, 48.0),
        ("C", "BB", "100-300", 150.0, 52.0),
        ("D", "BB", ">=1000", 1500.0, 58.0),
    ]
    multipliers = {
        2020: 1.00,
        2021: 1.10,
        2022: 1.02,
        2023: 1.03,
    }
    for nuts_id, country, band, density, latitude in specs:
        base = 1.0 + density / 1000.0
        for year in YEARS:
            mult = multipliers[year]
            rows.append(
                {
                    "NUTS_ID": nuts_id,
                    "CNTR_CODE": country,
                    "density_band": band,
                    "reference_population_density": density,
                    "centroid_latitude": latitude,
                    "year": year,
                    "radiance_mean_raw": base * mult,
                    "median_masked_mean": base * mult * 0.98,
                    "cf_cvg_mean": 40.0 + (year - 2020),
                }
            )
    synthetic = pd.DataFrame(rows)
    trans = make_transition_table(synthetic)
    annual = annual_summary(trans)
    density = density_summary(trans)
    diagnostics = diagnostic_table(trans)
    net = net_2020_2022(synthetic)

    if list(annual["end_year"]) != [2021, 2022, 2023]:
        raise RuntimeError("Self-test failed: unexpected annual transition years.")
    if len(diagnostics) != 3:
        raise RuntimeError("Self-test failed: diagnostics row count.")
    if len(net) != 4:
        raise RuntimeError("Self-test failed: net density row count.")
    if density.empty:
        raise RuntimeError("Self-test failed: density summary is empty.")
    required_net = {
        "net_2020_2022_geometric_pct_average_masked",
        "net_2020_2022_geometric_pct_median_masked",
    }
    if not required_net.issubset(net.columns):
        raise RuntimeError("Self-test failed: net endpoint columns missing.")
    print("Downstream synthetic self-test: PASS", flush=True)


def make_plots(
    trans: pd.DataFrame,
    annual: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    ax.plot(
        annual["end_year"],
        annual["average_masked_geometric_pct"],
        marker="o",
        label="average_masked",
    )
    ax.plot(
        annual["end_year"],
        annual["median_masked_geometric_pct"],
        marker="o",
        label="median_masked",
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Transition ending in year")
    ax.set_ylabel("Mean within-region radiance change (%)")
    ax.set_title("Annual VNL layer comparison on fixed NUTS2016 geography")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "average_vs_median_masked_2020_2023.pdf", dpi=dpi)
    plt.close(fig)

    g = trans[trans["year"] == 2021].copy()
    fig, ax = plt.subplots(figsize=(7.0, 6.4))
    ax.scatter(
        g["pct_change_average_masked"],
        g["pct_change_median_masked"],
        s=12,
        alpha=0.45,
    )
    lim = float(
        np.nanpercentile(
            np.abs(
                np.r_[
                    g["pct_change_average_masked"].values,
                    g["pct_change_median_masked"].values,
                ]
            ),
            99.0,
        )
    )
    lim = max(5.0, lim)
    ax.plot([-lim, lim], [-lim, lim], linestyle="--", linewidth=1.0)
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("average_masked change 2020→2021 (%)")
    ax.set_ylabel("median_masked change 2020→2021 (%)")
    ax.set_title("Does the 2021 pulse survive a robust annual summary?")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "average_vs_median_2021_region_change.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.3, 6.2))
    ax.scatter(
        g["delta_cf_cvg"],
        g["pct_change_average_masked"],
        s=12,
        alpha=0.45,
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("Change in mean cloud-free observation count, 2020→2021")
    ax.set_ylabel("average_masked radiance change (%)")
    ax.set_title("2021 radiance pulse versus cloud-free coverage")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "radiance_2021_vs_cf_cvg_change.pdf", dpi=dpi)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    config = load_config(args.config)
    diag_cfg = config["analysis"]["fixed_nuts2016_2021_diagnostic"]

    # Do not start eight expensive global-raster extractions until every
    # downstream pandas transformation has passed a synthetic smoke test.
    downstream_self_test()

    panel = add_latitude(load_locked_panel(config), config)
    manifest = locate_companions(panel, config)
    companions = extract_all(panel, manifest, config)
    analysis = build_analysis(panel, companions)
    trans = make_transition_table(analysis)

    annual = annual_summary(trans)
    density = density_summary(trans)
    diagnostics = diagnostic_table(trans)
    net = net_2020_2022(analysis)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "viirs_annual_product_layer_audit"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    outputs = {
        "manifest": processed / "viirs_annual_product_layer_manifest_2020_2023.csv",
        "regional": processed / "viirs_annual_product_layers_nuts2016_2020_2023.csv",
        "transitions": processed / "viirs_annual_product_layer_transitions_2021_2023.csv",
        "annual": processed / "viirs_annual_product_layer_annual_summary_2021_2023.csv",
        "density": processed / "viirs_annual_product_layer_density_summary_2021_2023.csv",
        "diagnostics": processed / "viirs_annual_product_layer_diagnostics_2021_2023.csv",
        "net": processed / "viirs_annual_product_layer_net_2020_2022_by_density.csv",
        "summary": processed / "viirs_annual_product_layer_audit_summary.txt",
    }

    manifest.to_csv(outputs["manifest"], index=False)
    analysis.to_csv(outputs["regional"], index=False)
    trans.to_csv(outputs["transitions"], index=False)
    annual.to_csv(outputs["annual"], index=False)
    density.to_csv(outputs["density"], index=False)
    diagnostics.to_csv(outputs["diagnostics"], index=False)
    net.to_csv(outputs["net"], index=False)

    make_plots(
        trans,
        annual,
        figures,
        int(diag_cfg["figure_dpi"]),
    )

    a21 = annual.loc[annual["end_year"] == 2021].iloc[0]
    d21 = diagnostics.loc[diagnostics["end_year"] == 2021].iloc[0]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("VIIRS ANNUAL PRODUCT-LAYER AUDIT: 2020-2023\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            "Purpose: determine whether the 2021 average_masked pulse is "
            "robust to median annual summarisation and independent of "
            "cloud-free observation coverage.\n"
        )
        h.write(
            f"Locked cohort: {analysis['NUTS_ID'].nunique()} NUTS3 regions; "
            f"{analysis['CNTR_CODE'].nunique()} countries.\n\n"
        )
        h.write("COMPANION RASTER MANIFEST\n")
        h.write(manifest.to_string(index=False))
        h.write("\n\nANNUAL LAYER COMPARISON\n")
        h.write(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nDENSITY-BAND COMPARISON\n")
        h.write(density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nREGION-LEVEL DIAGNOSTICS\n")
        h.write(
            diagnostics.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )
        h.write("\n\nNET 2020->2022 BY DENSITY BAND\n")
        h.write(net.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nKEY 2021 CHECKS\n")
        h.write(
            f"average_masked geometric change 2020->2021: "
            f"{a21['average_masked_geometric_pct']:+.4f}%\n"
        )
        h.write(
            f"median_masked geometric change 2020->2021: "
            f"{a21['median_masked_geometric_pct']:+.4f}%\n"
        )
        h.write(
            f"Spearman region-level average-vs-median change: "
            f"rho={d21['rho_average_vs_median_change']:.4f}\n"
        )
        h.write(
            f"Spearman average change vs delta cf_cvg: "
            f"rho={d21['rho_average_change_vs_delta_cf_cvg']:.4f}\n"
        )
        h.write(
            "\nDecision rule:\n"
            "- If the 2021 pulse is strongly attenuated in median_masked, "
            "average_masked is sensitive to within-year extremes/sampling and "
            "should not be accepted as the canonical annual trend layer "
            "without further justification.\n"
            "- If the pulse persists in median_masked but tracks cf_cvg, the "
            "annual observation/compositing process remains a likely driver.\n"
            "- If it persists in both layers and is weakly related to cf_cvg, "
            "the signal is not explained by these two product-quality checks; "
            "the next audit should then target the annual lit/background mask "
            "and the V2.1-to-V2.2 boundary explicitly.\n"
        )

    print("\nVIIRS ANNUAL PRODUCT-LAYER AUDIT")
    print(f"Locked cohort: {analysis['NUTS_ID'].nunique()} NUTS3")
    print("\nANNUAL LAYER COMPARISON")
    print(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nREGION-LEVEL DIAGNOSTICS")
    print(diagnostics.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nNET 2020 -> 2022 BY DENSITY BAND")
    print(net.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
