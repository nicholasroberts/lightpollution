#!/usr/bin/env python3
"""Validate the canonical density response against real historical VIIRS change.

This stage addresses a specific forecasting concern: the canonical future
response implies small negative annual changes in calibrated VIIRS DNB radiance
for dense NUTS3 regions. Before carrying that response into 2030-2050
projections, show directly whether the historical data support that pattern.

The validation deliberately uses one fixed geography, NUTS2016, over 2013-2020
because those years are all represented on the same NUTS release in the
production pipeline. It therefore avoids creating artificial regional slopes
across NUTS boundary changes.

Outputs include:
1. real, complete-panel historical trajectories by fixed population-density
   class, normalized to each NUTS3 region's 2013 radiance;
2. distributions of individual NUTS3 historical slopes by density class;
3. a Europe map of the observed NUTS3 slope;
4. the canonical density-response prediction for the same NUTS3 regions; and
5. the observed-minus-predicted residual map.

Important scale statement
-------------------------
The response is a NUTS3 mean calibrated VIIRS DNB radiance trend. It is NOT a
claim that the brightest pixels in city centres themselves dim at the same
rate. Dense NUTS3 regions include heterogeneous urban and peri-urban land.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from scipy.stats import linregress
from shapely.geometry import box

from project_config import configured_path, load_config
from production_validation import validate_production_merged


def load_fixed_panel(config: dict) -> pd.DataFrame:
    cfg = config["analysis"]["future_radiance_historical_validation"]
    processed = configured_path(config, "processed")
    path = processed / config["outputs"]["merged"]
    data = pd.read_parquet(path).copy()
    validate_production_merged(data, config)

    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    release = int(cfg["nuts_release"])
    years = list(range(first_year, last_year + 1))

    d = data[
        data["year"].isin(years)
        & (data["nuts_release"].astype(int) == release)
    ].copy()

    for col in ["population_density", "radiance_mean_corrected"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    good = (
        d["NUTS_ID"].notna()
        & d["CNTR_CODE"].notna()
        & np.isfinite(d["population_density"])
        & np.isfinite(d["radiance_mean_corrected"])
        & (d["population_density"] > 0)
        & (d["radiance_mean_corrected"] > 0)
    )
    d = d.loc[good].copy()

    duplicates = d.duplicated(["NUTS_ID", "year"], keep=False)
    if duplicates.any():
        raise ValueError(
            "Duplicate NUTS_ID-year rows in fixed historical validation panel."
        )

    counts = d.groupby("NUTS_ID")["year"].nunique()
    if bool(cfg["require_complete_panel"]):
        keep_ids = counts[counts == len(years)].index
    else:
        keep_ids = counts[
            counts >= int(cfg["minimum_panel_years"])
        ].index
    d = d[d["NUTS_ID"].isin(keep_ids)].copy()

    country_counts = d.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    keep_countries = country_counts[
        country_counts >= int(cfg["min_country_regions"])
    ].index
    d = d[d["CNTR_CODE"].isin(keep_countries)].copy()

    d["log10_radiance"] = np.log10(
        d["radiance_mean_corrected"].astype(float)
    )

    ref = (
        d.groupby("NUTS_ID", as_index=False)["population_density"]
        .median()
        .rename(
            columns={
                "population_density": "reference_population_density"
            }
        )
    )
    d = d.merge(ref, on="NUTS_ID", how="left", validate="many_to_one")

    density_cfg = config["analysis"]["density_temporal_change"]
    edges = (
        [-np.inf]
        + list(density_cfg["density_band_breaks"])
        + [np.inf]
    )
    d["density_band"] = pd.cut(
        d["reference_population_density"],
        bins=edges,
        labels=density_cfg["density_band_labels"],
        right=False,
        ordered=True,
    )

    expected = set(years)
    for nuts_id, g in d.groupby("NUTS_ID"):
        if bool(cfg["require_complete_panel"]) and set(g["year"]) != expected:
            raise RuntimeError(
                f"Incomplete panel survived filtering for {nuts_id}."
            )

    return d.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)


def region_slopes(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for nuts_id, g in panel.groupby("NUTS_ID", sort=False):
        g = g.sort_values("year")
        fit = linregress(g["year"], g["log10_radiance"])
        rows.append(
            {
                "NUTS_ID": nuts_id,
                "CNTR_CODE": str(g["CNTR_CODE"].iloc[0]),
                "density_band": str(g["density_band"].iloc[0]),
                "reference_population_density": float(
                    g["reference_population_density"].iloc[0]
                ),
                "n_years": int(len(g)),
                "observed_log10_slope_per_year": float(fit.slope),
                "observed_slope_se": float(fit.stderr),
                "observed_slope_p": float(fit.pvalue),
                "observed_slope_r2": float(fit.rvalue**2),
                "observed_pct_change_per_year": float(
                    100.0 * (10.0 ** fit.slope - 1.0)
                ),
            }
        )
    return pd.DataFrame(rows)


def add_canonical_predictions(
    slopes: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    processed = configured_path(config, "processed")
    response_path = processed / "future_radiance_canonical_response.csv"
    if not response_path.exists():
        raise FileNotFoundError(
            f"Missing {response_path}. Run future-radiance-response first."
        )

    response = pd.read_csv(response_path)
    required = [
        "population_density",
        "calibrated_viirs_change_pct_per_year",
        "annual_log10_calibrated_viirs_change",
    ]
    missing = [c for c in required if c not in response.columns]
    if missing:
        raise KeyError(
            f"Canonical response is missing required columns: {missing}"
        )

    if "selected_model" in response.columns:
        models = set(response["selected_model"].dropna().astype(str))
        if models != {"density_only"}:
            raise RuntimeError(
                "Historical validation currently expects the selected canonical "
                "response to be density_only."
            )

    response = response.sort_values("population_density").copy()
    x = np.log10(response["population_density"].to_numpy(dtype=float))
    y_log = response[
        "annual_log10_calibrated_viirs_change"
    ].to_numpy(dtype=float)

    log_d = np.log10(
        slopes["reference_population_density"].to_numpy(dtype=float)
    )
    pred_log = np.interp(log_d, x, y_log, left=y_log[0], right=y_log[-1])

    out = slopes.copy()
    out["predicted_log10_slope_per_year"] = pred_log
    out["predicted_pct_change_per_year"] = (
        100.0 * (10.0**pred_log - 1.0)
    )
    out["residual_pct_change_per_year"] = (
        out["observed_pct_change_per_year"]
        - out["predicted_pct_change_per_year"]
    )
    return out


def actual_class_trajectories(
    panel: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    first_year = int(
        config["analysis"]["future_radiance_historical_validation"][
            "first_year"
        ]
    )
    baseline = panel[panel["year"] == first_year][
        ["NUTS_ID", "log10_radiance"]
    ].rename(columns={"log10_radiance": "baseline_log10_radiance"})

    d = panel.merge(
        baseline,
        on="NUTS_ID",
        how="inner",
        validate="many_to_one",
    )
    d["delta_log10_from_baseline"] = (
        d["log10_radiance"] - d["baseline_log10_radiance"]
    )

    rows = []
    for (band, year), g in d.groupby(
        ["density_band", "year"], observed=True
    ):
        values = g["delta_log10_from_baseline"].to_numpy(dtype=float)
        mean_log = float(np.mean(values))
        se = (
            float(np.std(values, ddof=1) / np.sqrt(len(values)))
            if len(values) > 1
            else np.nan
        )
        rows.append(
            {
                "density_band": str(band),
                "year": int(year),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "mean_delta_log10_from_2013": mean_log,
                "se_delta_log10_from_2013": se,
                "observed_index_2013_100": 100.0 * 10.0**mean_log,
                "observed_index_ci_low": (
                    100.0 * 10.0 ** (mean_log - 1.96 * se)
                    if np.isfinite(se)
                    else np.nan
                ),
                "observed_index_ci_high": (
                    100.0 * 10.0 ** (mean_log + 1.96 * se)
                    if np.isfinite(se)
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def class_summary(
    slopes: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    processed = configured_path(config, "processed")
    headline = pd.read_csv(
        processed / "future_radiance_headline_density_classes.csv"
    )
    headline = headline[
        [
            "density_band",
            "reference_population_density_median",
            "pct_radiance_change_per_year",
            "pct_change_per_year_ci_low",
            "pct_change_per_year_ci_high",
        ]
    ].copy()
    headline["density_band"] = headline["density_band"].astype(str)

    rows = []
    labels = config["analysis"]["density_temporal_change"][
        "density_band_labels"
    ]
    for band in labels:
        g = slopes[slopes["density_band"] == str(band)].copy()
        if g.empty:
            continue
        obs = g["observed_pct_change_per_year"].to_numpy(dtype=float)
        pred = g["predicted_pct_change_per_year"].to_numpy(dtype=float)
        rows.append(
            {
                "density_band": str(band),
                "n_nuts3_fixed_2013_2020": int(len(g)),
                "fixed_panel_reference_density_median": float(
                    g["reference_population_density"].median()
                ),
                "observed_region_slope_mean_pct_per_year": float(
                    np.mean(obs)
                ),
                "observed_region_slope_median_pct_per_year": float(
                    np.median(obs)
                ),
                "observed_region_slope_q25": float(np.quantile(obs, 0.25)),
                "observed_region_slope_q75": float(np.quantile(obs, 0.75)),
                "share_regions_negative_slope": float(np.mean(obs < 0)),
                "canonical_prediction_mean_pct_per_year": float(
                    np.mean(pred)
                ),
                "mean_observed_minus_predicted_pct_per_year": float(
                    np.mean(obs - pred)
                ),
            }
        )

    out = pd.DataFrame(rows)
    out = out.merge(
        headline,
        on="density_band",
        how="left",
        validate="one_to_one",
    )
    return out


def plot_class_trajectories(
    trajectories: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict,
    output: Path,
) -> None:
    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    first_year = int(
        config["analysis"]["future_radiance_historical_validation"][
            "first_year"
        ]
    )

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), sharex=True)
    axes = axes.ravel()

    for i, band in enumerate(labels):
        ax = axes[i]
        g = trajectories[
            trajectories["density_band"] == str(band)
        ].sort_values("year")
        row = summary[summary["density_band"] == str(band)].iloc[0]

        years = g["year"].to_numpy(dtype=float)
        ax.plot(
            years,
            g["observed_index_2013_100"],
            marker="o",
            linewidth=1.5,
            label="Observed complete-panel data",
        )
        ax.fill_between(
            years,
            g["observed_index_ci_low"],
            g["observed_index_ci_high"],
            alpha=0.18,
        )

        canonical_rate = float(row["pct_radiance_change_per_year"])
        model_index = 100.0 * (
            1.0 + canonical_rate / 100.0
        ) ** (years - first_year)
        ax.plot(
            years,
            model_index,
            linestyle="--",
            linewidth=1.8,
            label="Canonical class trend",
        )

        ax.axhline(100, linewidth=0.7)
        ax.set_title(
            f"{band} people km$^{{-2}}$\n"
            f"canonical {canonical_rate:+.2f}% yr$^{{-1}}$"
        )
        ax.grid(alpha=0.2)

    axes[5].axis("off")
    handles, legend_labels = axes[0].get_legend_handles_labels()
    axes[5].legend(
        handles,
        legend_labels,
        loc="center",
        frameon=False,
        fontsize=10,
    )

    fig.supylabel("Calibrated VIIRS DNB radiance index (2013 = 100)")
    fig.supxlabel("Year")
    fig.suptitle(
        "Real historical radiance trajectories versus canonical density-class trends\n"
        "Fixed NUTS2016 complete panel, 2013–2020",
        fontsize=14,
    )
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.94))
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_slope_distributions(
    slopes: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict,
    output: Path,
) -> None:
    labels = list(
        config["analysis"]["density_temporal_change"][
            "density_band_labels"
        ]
    )
    values = [
        slopes.loc[
            slopes["density_band"] == str(b),
            "observed_pct_change_per_year",
        ].to_numpy(dtype=float)
        for b in labels
    ]

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.boxplot(
        values,
        tick_labels=labels,
        showfliers=False,
        widths=0.62,
        patch_artist=False,
    )

    for i, band in enumerate(labels, start=1):
        row = summary[summary["density_band"] == str(band)].iloc[0]
        ax.scatter(
            i,
            row["pct_radiance_change_per_year"],
            marker="D",
            s=52,
            zorder=4,
            label=(
                "Canonical density-class estimate"
                if i == 1
                else None
            ),
        )

    ax.axhline(0, linewidth=1)
    ax.set_xlabel("Fixed population-density class (people km$^{-2}$)")
    ax.set_ylabel(
        "Observed NUTS3 calibrated VIIRS DNB trend (% yr$^{-1}$)"
    )
    ax.set_title(
        "Distribution of real 2013–2020 NUTS3 radiance trends\n"
        "diamond = canonical 2013–2024 class estimate"
    )
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def load_geometry(config: dict) -> gpd.GeoDataFrame:
    cfg = config["analysis"]["future_radiance_historical_validation"]
    release = int(cfg["nuts_release"])
    nuts_raw = configured_path(config, "nuts_raw")
    path = nuts_raw / config["nuts"]["filename_template"].format(
        release=release
    )
    if not path.exists():
        raise FileNotFoundError(path)

    geometry = gpd.read_file(path)
    west, south, east, north = map(
        float, config["analysis"]["europe_maps"]["extent_4326"]
    )
    europe_box = gpd.GeoDataFrame(
        geometry=[box(west, south, east, north)],
        crs="EPSG:4326",
    )
    geometry = gpd.clip(geometry, europe_box)
    return geometry, europe_box


def plot_maps(
    slopes: pd.DataFrame,
    config: dict,
    output: Path,
) -> None:
    geometry, europe_box = load_geometry(config)
    id_field = config["nuts"]["id_field"]
    mapped = geometry.merge(
        slopes[
            [
                "NUTS_ID",
                "observed_pct_change_per_year",
                "predicted_pct_change_per_year",
                "residual_pct_change_per_year",
                "density_band",
            ]
        ],
        left_on=id_field,
        right_on="NUTS_ID",
        how="left",
        validate="one_to_one",
    ).to_crs("EPSG:3035")

    box3035 = europe_box.to_crs("EPSG:3035")
    xmin, ymin, xmax, ymax = box3035.total_bounds

    cfg = config["analysis"]["future_radiance_historical_validation"]
    q = float(cfg["map_absolute_quantile"])

    common_values = np.concatenate(
        [
            slopes["observed_pct_change_per_year"].to_numpy(dtype=float),
            slopes["predicted_pct_change_per_year"].to_numpy(dtype=float),
        ]
    )
    common_values = common_values[np.isfinite(common_values)]
    common_lim = float(np.quantile(np.abs(common_values), q))
    common_lim = max(common_lim, 0.25)

    residual_values = slopes[
        "residual_pct_change_per_year"
    ].to_numpy(dtype=float)
    residual_values = residual_values[np.isfinite(residual_values)]
    residual_lim = float(np.quantile(np.abs(residual_values), q))
    residual_lim = max(residual_lim, 0.25)

    cmap = "coolwarm"
    common_norm = TwoSlopeNorm(
        vmin=-common_lim, vcenter=0, vmax=common_lim
    )
    residual_norm = TwoSlopeNorm(
        vmin=-residual_lim, vcenter=0, vmax=residual_lim
    )

    panels = [
        (
            "observed_pct_change_per_year",
            "Observed historical trend\n2013–2020 fixed NUTS2016",
            common_norm,
        ),
        (
            "predicted_pct_change_per_year",
            "Canonical density-response prediction\nfor the same NUTS3 regions",
            common_norm,
        ),
        (
            "residual_pct_change_per_year",
            "Residual: observed − predicted",
            residual_norm,
        ),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17.0, 7.0))
    for ax, (column, title, norm) in zip(axes, panels):
        mapped.plot(
            ax=ax,
            color="0.92",
            edgecolor="0.78",
            linewidth=0.15,
        )
        good = mapped[np.isfinite(mapped[column])].copy()
        good.plot(
            column=column,
            ax=ax,
            cmap=cmap,
            norm=norm,
            edgecolor="white",
            linewidth=0.08,
        )
        mapped.dissolve(by="CNTR_CODE").boundary.plot(
            ax=ax,
            color="0.18",
            linewidth=0.35,
        )
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_axis_off()
        ax.set_title(title, fontsize=12)

    sm1 = ScalarMappable(norm=common_norm, cmap=cmap)
    sm1.set_array([])
    cbar1 = fig.colorbar(
        sm1,
        ax=axes[:2].tolist(),
        orientation="horizontal",
        fraction=0.045,
        pad=0.035,
    )
    cbar1.set_label(
        "Annual change in calibrated VIIRS DNB radiance (% yr$^{-1}$)"
    )

    sm2 = ScalarMappable(norm=residual_norm, cmap=cmap)
    sm2.set_array([])
    cbar2 = fig.colorbar(
        sm2,
        ax=[axes[2]],
        orientation="horizontal",
        fraction=0.045,
        pad=0.035,
    )
    cbar2.set_label("Observed − predicted (% yr$^{-1}$)")

    fig.suptitle(
        "Does the canonical population-density response reproduce the real "
        "historical spatial pattern?",
        fontsize=15,
        y=0.97,
    )
    fig.text(
        0.5,
        0.015,
        "Observed slopes are region-wise OLS fits to annual production-calibrated "
        "VIIRS DNB radiance. Fixed NUTS2016 geography avoids boundary-change artefacts.",
        ha="center",
        fontsize=9,
    )
    fig.savefig(output, dpi=int(cfg["figure_dpi"]), bbox_inches="tight")
    plt.close(fig)


def write_summary(
    path: Path,
    panel: pd.DataFrame,
    slopes: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict,
) -> None:
    cfg = config["analysis"]["future_radiance_historical_validation"]

    obs = slopes["observed_pct_change_per_year"].to_numpy(dtype=float)
    pred = slopes["predicted_pct_change_per_year"].to_numpy(dtype=float)
    finite = np.isfinite(obs) & np.isfinite(pred)

    corr = (
        float(np.corrcoef(obs[finite], pred[finite])[0, 1])
        if finite.sum() > 2
        else np.nan
    )
    mae = float(np.mean(np.abs(obs[finite] - pred[finite])))
    bias = float(np.mean(obs[finite] - pred[finite]))

    with path.open("w", encoding="utf-8") as h:
        h.write(
            "HISTORICAL SPATIAL VALIDATION OF THE CANONICAL VIIRS DENSITY RESPONSE\n"
        )
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Fixed geography: NUTS {cfg['nuts_release']}; "
            f"years {cfg['first_year']}-{cfg['last_year']}.\n"
        )
        h.write(
            f"Complete-panel NUTS3: {slopes['NUTS_ID'].nunique()}; "
            f"countries: {slopes['CNTR_CODE'].nunique()}; "
            f"annual observations: {len(panel)}.\n"
        )
        h.write(
            "Observed NUTS3 trend = OLS slope of log10 production-calibrated "
            "VIIRS DNB radiance against year, converted to % yr-1.\n"
        )
        h.write(
            "Predicted NUTS3 trend = canonical density-only response evaluated "
            "at that region's median 2013-2020 population density.\n\n"
        )

        h.write("DENSITY-CLASS VALIDATION\n")
        h.write(
            summary.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nOVERALL REGION-LEVEL AGREEMENT\n")
        h.write(
            f"Observed-predicted Pearson correlation: {corr:.6g}\n"
            f"Mean absolute error: {mae:.6g} percentage points yr-1\n"
            f"Mean observed-minus-predicted bias: {bias:.6g} "
            "percentage points yr-1\n"
        )
        h.write("\nINTERPRETATION LIMIT\n")
        h.write(
            "The map validates a NUTS3-mean regional trend. It does not imply "
            "that city-centre pixels themselves decline at the regional mean "
            "rate. The >=1000 people km-2 class comprises dense NUTS3 regions "
            "with heterogeneous urban and peri-urban land.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    panel = load_fixed_panel(config)
    slopes = region_slopes(panel)
    slopes = add_canonical_predictions(slopes, config)
    trajectories = actual_class_trajectories(panel, config)
    summary = class_summary(slopes, config)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "future_radiance_historical_validation"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    paths = {
        "slopes": processed
        / "future_radiance_historical_validation_nuts3.csv",
        "trajectories": processed
        / "future_radiance_historical_validation_trajectories.csv",
        "summary_table": processed
        / "future_radiance_historical_validation_by_density.csv",
        "summary": processed
        / "future_radiance_historical_validation_summary.txt",
    }

    slopes.to_csv(paths["slopes"], index=False)
    trajectories.to_csv(paths["trajectories"], index=False)
    summary.to_csv(paths["summary_table"], index=False)

    fig_trajectories = figures / "historical_real_vs_canonical_trajectories.pdf"
    fig_distribution = figures / "historical_nuts3_slope_distributions.pdf"
    fig_maps = figures / "historical_observed_predicted_residual_maps.pdf"

    plot_class_trajectories(
        trajectories, summary, config, fig_trajectories
    )
    plot_slope_distributions(
        slopes, summary, config, fig_distribution
    )
    plot_maps(slopes, config, fig_maps)
    write_summary(
        paths["summary"], panel, slopes, summary, config
    )

    print("\nHISTORICAL DENSITY-CLASS VALIDATION")
    print(
        summary.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nWROTE")
    for p in [
        *paths.values(),
        fig_trajectories,
        fig_distribution,
        fig_maps,
    ]:
        print(p)


if __name__ == "__main__":
    main()
