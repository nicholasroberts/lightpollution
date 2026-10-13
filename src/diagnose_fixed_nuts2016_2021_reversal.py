#!/usr/bin/env python3
"""Diagnose the 2020->2021 pulse and 2021->2022 reversal in fixed NUTS2016 VIIRS.

This stage is deliberately diagnostic. It reuses the locked fixed-NUTS2016
2014-2024 extraction and does not alter the production calibration, cohort,
density bands, or any upstream files.

Questions tested
----------------
1. Is the 2021 pulse already present in the raw annual radiance, or created by
   the production 0.235 nW cm-2 sr-1 correction?
2. Does the 2022 change reverse the 2021 change within the same NUTS3 regions?
3. Is the pulse/reversal associated with latitude and/or baseline brightness?
4. After the two-year 2020->2022 interval, how much net change remains by
   density band and country?
5. Which exact local annual raster file supplied each year?

The point is to distinguish a transient annual-product/sampling discontinuity
from a persistent ecological/anthropogenic trend before fitting any final
historical model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr

from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import (
    load_fixed_nuts,
    prepare_analysis_panel,
)


TRANSITION_YEARS = [2020, 2021, 2022, 2023]


def _restore_density_order(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    panel = panel.copy()
    panel["density_band"] = pd.Categorical(
        panel["density_band"],
        categories=labels,
        ordered=True,
    )
    if panel["density_band"].isna().any():
        raise RuntimeError(
            "Density-band labels do not match configured ordering: "
            + repr(labels)
        )
    return panel


def load_locked_panel(config: dict) -> pd.DataFrame:
    cfg = config["analysis"]["fixed_nuts2016_2014_2024"]
    diag = config["analysis"]["fixed_nuts2016_2021_diagnostic"]
    processed = configured_path(config, "processed")

    fixed_path = processed / diag["fixed_radiance_file"]
    reference_path = processed / diag["population_reference_file"]
    if not fixed_path.exists() or not reference_path.exists():
        raise FileNotFoundError(
            "Missing locked fixed-NUTS2016 outputs. Run "
            "'python src/run_pipeline.py future-radiance-fixed-2014-2024' first."
        )

    fixed = pd.read_parquet(fixed_path).copy()
    reference = pd.read_csv(reference_path).copy()
    panel = prepare_analysis_panel(fixed, reference, config)
    panel = _restore_density_order(panel, config)
    panel = panel.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)

    expected_n = int(diag["expected_n_nuts3"])
    if panel["NUTS_ID"].nunique() != expected_n:
        raise RuntimeError(
            f"Expected {expected_n} fixed NUTS3 regions, found "
            f"{panel['NUTS_ID'].nunique()}."
        )

    expected_years = list(range(int(cfg["first_year"]), int(cfg["last_year"]) + 1))
    observed_years = sorted(panel["year"].astype(int).unique().tolist())
    if observed_years != expected_years:
        raise RuntimeError(
            f"Expected years {expected_years}; observed {observed_years}."
        )

    counts = panel.groupby("NUTS_ID")["year"].nunique()
    if not (counts == len(expected_years)).all():
        raise RuntimeError("Locked panel is no longer complete by NUTS_ID.")

    return panel


def add_latitude(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release)[["NUTS_ID", "geometry"]].copy()

    # Centroids are computed in the European equal-area CRS, then transformed
    # back to WGS84 to obtain a meaningful representative latitude.
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
            "centroid_longitude": centres.geometry.x.to_numpy(float),
        }
    )
    return panel.merge(loc, on="NUTS_ID", how="left", validate="many_to_one")


def make_differences(panel: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "NUTS_ID",
        "CNTR_CODE",
        "density_band",
        "reference_population_density",
        "centroid_latitude",
        "centroid_longitude",
        "year",
        "radiance_mean_raw",
        "radiance_mean_corrected",
        "viirs_file",
    ]
    d = panel[cols].copy().sort_values(["NUTS_ID", "year"])

    for metric in ["raw", "corrected"]:
        col = f"radiance_mean_{metric}"
        d[col] = pd.to_numeric(d[col], errors="coerce")
        if (~np.isfinite(d[col]) | (d[col] <= 0)).any():
            raise RuntimeError(
                f"Non-positive or non-finite {col} in locked analysis panel."
            )
        d[f"log10_{metric}"] = np.log10(d[col])
        d[f"previous_log10_{metric}"] = d.groupby("NUTS_ID")[
            f"log10_{metric}"
        ].shift(1)
        d[f"delta_log10_{metric}"] = (
            d[f"log10_{metric}"] - d[f"previous_log10_{metric}"]
        )
        d[f"pct_change_{metric}"] = 100.0 * (
            10.0 ** d[f"delta_log10_{metric}"] - 1.0
        )

    d["previous_year"] = d.groupby("NUTS_ID")["year"].shift(1)
    d["previous_raw"] = d.groupby("NUTS_ID")["radiance_mean_raw"].shift(1)
    d["previous_corrected"] = d.groupby("NUTS_ID")[
        "radiance_mean_corrected"
    ].shift(1)
    d = d[d["previous_year"].notna()].copy()
    if (
        d["year"].astype(int) - d["previous_year"].astype(int) != 1
    ).any():
        raise RuntimeError("Non-consecutive years found in locked panel.")

    d["end_year"] = d["year"].astype(int)
    return d


def geom_pct(delta_log10: pd.Series) -> float:
    return float(100.0 * (10.0 ** delta_log10.mean() - 1.0))


def transition_summary(diff: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for year, g in diff.groupby("end_year", observed=True):
        if int(year) not in TRANSITION_YEARS:
            continue
        row = {
            "end_year": int(year),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
        }
        for metric in ["raw", "corrected"]:
            x = g[f"delta_log10_{metric}"]
            row[f"geometric_pct_change_{metric}"] = geom_pct(x)
            row[f"median_region_pct_change_{metric}"] = float(
                g[f"pct_change_{metric}"].median()
            )
            row[f"pct_regions_increasing_{metric}"] = float(
                100.0 * (x > 0).mean()
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values("end_year")


def density_transition_summary(diff: pd.DataFrame) -> pd.DataFrame:
    rows = []
    chosen = diff[diff["end_year"].isin([2021, 2022])].copy()
    for (year, band), g in chosen.groupby(
        ["end_year", "density_band"], observed=True
    ):
        row = {
            "end_year": int(year),
            "density_band": str(band),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
        }
        for metric in ["raw", "corrected"]:
            row[f"geometric_pct_change_{metric}"] = geom_pct(
                g[f"delta_log10_{metric}"]
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["end_year", "density_band"])


def make_reversal_table(panel: pd.DataFrame) -> pd.DataFrame:
    keep = panel[panel["year"].isin([2020, 2021, 2022])].copy()
    id_cols = [
        "NUTS_ID",
        "CNTR_CODE",
        "density_band",
        "reference_population_density",
        "centroid_latitude",
        "centroid_longitude",
    ]
    base = keep[id_cols].drop_duplicates("NUTS_ID")
    if base["NUTS_ID"].nunique() != panel["NUTS_ID"].nunique():
        raise RuntimeError("Failed to retain all fixed NUTS3 regions.")

    wide = base.copy()
    for metric in ["raw", "corrected"]:
        p = keep.pivot(
            index="NUTS_ID",
            columns="year",
            values=f"radiance_mean_{metric}",
        )
        for year in [2020, 2021, 2022]:
            wide = wide.merge(
                p[[year]].rename(columns={year: f"{metric}_{year}"}),
                left_on="NUTS_ID",
                right_index=True,
                how="left",
                validate="one_to_one",
            )

        wide[f"d21_{metric}"] = np.log10(
            wide[f"{metric}_2021"] / wide[f"{metric}_2020"]
        )
        wide[f"d22_{metric}"] = np.log10(
            wide[f"{metric}_2022"] / wide[f"{metric}_2021"]
        )
        wide[f"net20_22_{metric}"] = np.log10(
            wide[f"{metric}_2022"] / wide[f"{metric}_2020"]
        )
        wide[f"pct21_{metric}"] = 100.0 * (
            10.0 ** wide[f"d21_{metric}"] - 1.0
        )
        wide[f"pct22_{metric}"] = 100.0 * (
            10.0 ** wide[f"d22_{metric}"] - 1.0
        )
        wide[f"net_pct20_22_{metric}"] = 100.0 * (
            10.0 ** wide[f"net20_22_{metric}"] - 1.0
        )
    return wide


def _spearman(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    ok = np.isfinite(x.to_numpy(float)) & np.isfinite(y.to_numpy(float))
    r = spearmanr(x[ok], y[ok])
    return float(r.statistic), float(r.pvalue)


def reversal_diagnostics(rev: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in ["raw", "corrected"]:
        d21 = rev[f"d21_{metric}"]
        d22 = rev[f"d22_{metric}"]
        rho_reversal, p_reversal = _spearman(d21, d22)
        rho_lat21, p_lat21 = _spearman(rev["centroid_latitude"], d21)
        rho_lat22, p_lat22 = _spearman(rev["centroid_latitude"], d22)
        rho_bright21, p_bright21 = _spearman(
            np.log10(rev[f"{metric}_2020"]), d21
        )
        fit = linregress(d21.to_numpy(float), d22.to_numpy(float))
        rows.append(
            {
                "metric": metric,
                "n_nuts3": int(len(rev)),
                "spearman_d21_vs_d22": rho_reversal,
                "spearman_d21_vs_d22_p": p_reversal,
                "ols_d22_on_d21_slope": float(fit.slope),
                "ols_d22_on_d21_r2": float(fit.rvalue ** 2),
                "pct_opposite_sign_2021_2022": float(
                    100.0 * ((d21 * d22) < 0).mean()
                ),
                "pct_2021_up_2022_down": float(
                    100.0 * ((d21 > 0) & (d22 < 0)).mean()
                ),
                "spearman_latitude_vs_d21": rho_lat21,
                "latitude_vs_d21_p": p_lat21,
                "spearman_latitude_vs_d22": rho_lat22,
                "latitude_vs_d22_p": p_lat22,
                "spearman_log2020_brightness_vs_d21": rho_bright21,
                "brightness_vs_d21_p": p_bright21,
                "geometric_2020_2022_net_pct": float(
                    100.0 * (10.0 ** rev[f"net20_22_{metric}"].mean() - 1.0)
                ),
            }
        )
    return pd.DataFrame(rows)


def net_by_density(rev: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band, g in rev.groupby("density_band", observed=True):
        row = {
            "density_band": str(band),
            "n_nuts3": int(len(g)),
        }
        for metric in ["raw", "corrected"]:
            row[f"net_2020_2022_geometric_pct_{metric}"] = float(
                100.0
                * (
                    10.0 ** g[f"net20_22_{metric}"].mean()
                    - 1.0
                )
            )
            row[f"median_region_net_pct_{metric}"] = float(
                g[f"net_pct20_22_{metric}"].median()
            )
        rows.append(row)
    return pd.DataFrame(rows)


def country_reversal(rev: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for country, g in rev.groupby("CNTR_CODE", observed=True):
        row = {
            "CNTR_CODE": str(country),
            "n_nuts3": int(len(g)),
            "median_latitude": float(g["centroid_latitude"].median()),
        }
        for metric in ["raw", "corrected"]:
            row[f"median_2021_pct_{metric}"] = float(
                g[f"pct21_{metric}"].median()
            )
            row[f"median_2022_pct_{metric}"] = float(
                g[f"pct22_{metric}"].median()
            )
            row[f"median_net_2020_2022_pct_{metric}"] = float(
                g[f"net_pct20_22_{metric}"].median()
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values("median_latitude")


def source_manifest(panel: pd.DataFrame, config: dict) -> pd.DataFrame:
    rows = []
    raw_dir = configured_path(config, "viirs_raw")
    for year, g in panel.groupby("year", observed=True):
        files = sorted(set(g["viirs_file"].astype(str)))
        if len(files) != 1:
            raise RuntimeError(
                f"Expected one source VIIRS filename for {year}; found {files}"
            )
        filename = files[0]
        path = raw_dir / filename
        rows.append(
            {
                "year": int(year),
                "viirs_file": filename,
                "local_file_exists_at_raw_root": bool(path.exists()),
                "local_size_bytes": int(path.stat().st_size) if path.exists() else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values("year")


def raw_directory_inventory(config: dict) -> pd.DataFrame:
    root = configured_path(config, "viirs_raw")
    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if any(str(y) in path.name for y in [2020, 2021, 2022, 2023]):
            rows.append(
                {
                    "relative_path": str(path.relative_to(root)),
                    "size_bytes": int(path.stat().st_size),
                    "suffix": path.suffix.lower(),
                }
            )
    return pd.DataFrame(rows)


def make_plots(
    rev: pd.DataFrame,
    net_density: pd.DataFrame,
    labels: list[str],
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # 2021 pulse versus 2022 reversal, using raw radiance so the plot cannot
    # be an artefact of the additive production correction.
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    ax.scatter(rev["pct21_raw"], rev["pct22_raw"], s=12, alpha=0.45)
    lim = np.nanpercentile(
        np.abs(np.r_[rev["pct21_raw"].values, rev["pct22_raw"].values]),
        99.0,
    )
    lim = max(5.0, float(lim))
    ax.plot([-lim, lim], [lim, -lim], linestyle="--", linewidth=1.0)
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Raw radiance change 2020→2021 (%)")
    ax.set_ylabel("Raw radiance change 2021→2022 (%)")
    ax.set_title("Within-region 2021 pulse versus 2022 reversal")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "raw_2021_pulse_vs_2022_reversal.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    ax.scatter(
        rev["centroid_latitude"],
        rev["pct21_raw"],
        s=12,
        alpha=0.45,
        label="2020→2021",
    )
    ax.scatter(
        rev["centroid_latitude"],
        rev["pct22_raw"],
        s=12,
        alpha=0.45,
        label="2021→2022",
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("NUTS3 centroid latitude (°N)")
    ax.set_ylabel("Raw radiance change (%)")
    ax.set_title("Pulse/reversal versus latitude")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "raw_pulse_reversal_vs_latitude.pdf", dpi=dpi)
    plt.close(fig)

    lookup = net_density.set_index("density_band")
    vals = [
        lookup.loc[label, "net_2020_2022_geometric_pct_raw"]
        for label in labels
    ]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    x = np.arange(len(labels))
    ax.bar(x, vals)
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Fixed population-density band (people km$^{-2}$)")
    ax.set_ylabel("Raw net radiance change, 2020→2022 (%)")
    ax.set_title("Two-year net change after the 2021 pulse")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(output_dir / "raw_net_2020_2022_by_density_band.pdf", dpi=dpi)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    config = load_config(args.config)
    diag_cfg = config["analysis"]["fixed_nuts2016_2021_diagnostic"]
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )

    panel = load_locked_panel(config)
    panel = add_latitude(panel, config)
    diff = make_differences(panel)

    annual = transition_summary(diff)
    density = density_transition_summary(diff)
    rev = make_reversal_table(panel)
    diagnostics = reversal_diagnostics(rev)
    net_density = net_by_density(rev)
    countries = country_reversal(rev)
    manifest = source_manifest(panel, config)
    inventory = raw_directory_inventory(config)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_2021_reversal"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    outputs = {
        "annual": processed / "fixed_nuts2016_2020_2023_raw_corrected_changes.csv",
        "density": processed / "fixed_nuts2016_2021_2022_raw_corrected_by_density.csv",
        "reversal": processed / "fixed_nuts2016_2021_2022_region_reversal.csv",
        "diagnostics": processed / "fixed_nuts2016_2021_2022_reversal_diagnostics.csv",
        "net_density": processed / "fixed_nuts2016_2020_2022_net_by_density.csv",
        "countries": processed / "fixed_nuts2016_2021_2022_country_reversal.csv",
        "manifest": processed / "fixed_nuts2016_viirs_source_manifest_2014_2024.csv",
        "inventory": processed / "viirs_raw_inventory_2020_2023.csv",
        "summary": processed / "fixed_nuts2016_2021_2022_reversal_summary.txt",
    }

    annual.to_csv(outputs["annual"], index=False)
    density.to_csv(outputs["density"], index=False)
    rev.to_csv(outputs["reversal"], index=False)
    diagnostics.to_csv(outputs["diagnostics"], index=False)
    net_density.to_csv(outputs["net_density"], index=False)
    countries.to_csv(outputs["countries"], index=False)
    manifest.to_csv(outputs["manifest"], index=False)
    inventory.to_csv(outputs["inventory"], index=False)

    make_plots(
        rev,
        net_density,
        labels,
        figures,
        int(diag_cfg["figure_dpi"]),
    )

    raw_diag = diagnostics.loc[diagnostics["metric"] == "raw"].iloc[0]
    cor_diag = diagnostics.loc[diagnostics["metric"] == "corrected"].iloc[0]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 DIAGNOSTIC: 2021 PULSE AND 2022 REVERSAL\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Locked panel: {panel['NUTS_ID'].nunique()} NUTS3 regions; "
            f"{panel['CNTR_CODE'].nunique()} countries; "
            f"{panel['year'].min()}-{panel['year'].max()}.\n"
        )
        h.write(
            "No upstream data, calibration, cohort, or density-band definition "
            "was changed.\n\n"
        )
        h.write("SOURCE RASTER MANIFEST\n")
        h.write(manifest.to_string(index=False))
        h.write("\n\nRAW-DIRECTORY INVENTORY FOR 2020-2023\n")
        h.write(
            inventory.to_string(index=False)
            if len(inventory)
            else "(No matching files found in raw VIIRS directory.)"
        )
        h.write("\n\nANNUAL RAW VERSUS CORRECTED CHANGE\n")
        h.write(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\n2021/2022 CHANGE BY DENSITY BAND\n")
        h.write(density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nREVERSAL DIAGNOSTICS\n")
        h.write(
            diagnostics.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )
        h.write("\n\nNET 2020->2022 CHANGE BY DENSITY BAND\n")
        h.write(
            net_density.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )
        h.write("\n\nCOUNTRY MEDIAN PULSE/REVERSAL\n")
        h.write(
            countries.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )
        h.write("\n\nPRE-SPECIFIED INTERPRETATION CHECKS\n")
        h.write(
            f"- Raw 2021/2022 reversal Spearman rho: "
            f"{raw_diag['spearman_d21_vs_d22']:.4f}.\n"
        )
        h.write(
            f"- Raw fraction 2021-up then 2022-down: "
            f"{raw_diag['pct_2021_up_2022_down']:.1f}%.\n"
        )
        h.write(
            f"- Raw latitude association with 2021 change: "
            f"rho={raw_diag['spearman_latitude_vs_d21']:.4f}.\n"
        )
        h.write(
            f"- Raw 2020-brightness association with 2021 change: "
            f"rho={raw_diag['spearman_log2020_brightness_vs_d21']:.4f}.\n"
        )
        h.write(
            f"- Raw overall net 2020->2022 change: "
            f"{raw_diag['geometric_2020_2022_net_pct']:+.4f}%.\n"
        )
        h.write(
            f"- Corrected overall net 2020->2022 change: "
            f"{cor_diag['geometric_2020_2022_net_pct']:+.4f}%.\n"
        )
        h.write(
            "\nA strong negative 2021-versus-2022 association, especially if "
            "already present in RAW radiance and structured by latitude or "
            "baseline brightness, is evidence against treating 2021 as a "
            "persistent biological/anthropogenic trend. It would justify a "
            "separate product-layer check using cloud-free coverage and "
            "median-masked radiance before final historical modelling.\n"
        )

    print("\nLOCKED PANEL")
    print(
        f"{panel['NUTS_ID'].nunique()} NUTS3; "
        f"{panel['CNTR_CODE'].nunique()} countries; "
        f"{len(panel)} annual rows"
    )
    print("\nSOURCE RASTER MANIFEST")
    print(manifest.to_string(index=False))
    print("\nANNUAL RAW VERSUS CORRECTED CHANGE")
    print(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nREVERSAL DIAGNOSTICS")
    print(diagnostics.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nNET 2020 -> 2022 CHANGE BY DENSITY BAND")
    print(net_density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
