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
    if mean_col not in result.columns:
        raise RuntimeError(
            f"Expected {mean_col} from exactextract; found {list(result.columns)}"
        )
    result = result[result[id_field].isin(keep_ids)][[id_field, mean_col]].copy()
    result.insert(0, "year", year)
    return result


def extract_all(
    panel: pd.DataFrame,
    manifest: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release)
    keep_ids = set(panel["NUTS_ID"].astype(str))
    paths = companion_paths(manifest, config)

    pieces = []
    for year in YEARS:
        year_parts = []
        for layer in LAYERS:
            year_parts.append(
                extract_layer(
                    paths[(year, layer)],
                    layer,
                    year,
                    nuts,
                    keep_ids,
                    config,
                )
            )
        merged = year_parts[0]
        for part in year_parts[1:]:
            merged = merged.merge(
                part, on=["year", "NUTS_ID"], how="inner", validate="one_to_one"
            )
        pieces.append(merged)
    return pd.concat(pieces, ignore_index=True)


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
    x["pct_change_cf_cvg"] = 100.0 * (
        x["cf_cvg_mean"] / x["prev_cf_cvg"] - 1.0
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
    wide = d[d["year"].isin([2020, 2022])].pivot(
        index="NUTS_ID",
        columns="year",
        values=["radiance_mean_raw", "median_masked_mean"],
    )
    meta = d[
        ["NUTS_ID", "density_band", "CNTR_CODE"]
    ].drop_duplicates("NUTS_ID").set_index("NUTS_ID")
    wide = meta.join(wide)
    rows = []
    for band, g in wide.groupby("density_band", observed=True):
        row = {"density_band": str(band), "n_nuts3": int(len(g))}
        for label, metric in [
            ("average_masked", "radiance_mean_raw"),
            ("median_masked", "median_masked_mean"),
        ]:
            a = g[(metric, 2020)]
            b = g[(metric, 2022)]
            ok = (a > 0) & (b > 0) & np.isfinite(a) & np.isfinite(b)
            delta = np.log10(b[ok] / a[ok])
            row[f"net_2020_2022_geometric_pct_{label}"] = float(
                100.0 * (10.0 ** delta.mean() - 1.0)
            )
            row[f"n_valid_{label}"] = int(ok.sum())
        rows.append(row)
    return pd.DataFrame(rows)


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
