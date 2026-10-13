#!/usr/bin/env python3
"""Audit absolute year-to-year VIIRS offsets on fixed NUTS2016 geography.

Scientific purpose
------------------
The 2021 percentage/log-radiance pulse is much larger in dark, low-density
regions. A small additive radiance disturbance can generate exactly that
pattern because the same absolute increment is a much larger percentage of a
dim baseline. This stage therefore asks whether 2020->2021 is unusual in
ABSOLUTE radiance units, and whether its absolute change is approximately flat
across independently defined baseline-brightness strata.

Primary response:
    raw EOG average_masked NUTS3 radiance, nW cm-2 sr-1

Sensitivity/positive control:
    production-corrected radiance. The known 2017 calibration discontinuity
    should be attenuated by the production correction; 2021 should not be
    created by that correction because it is applied on both sides of the
    2020->2021 transition.

Baseline brightness is fixed once per region as the median RAW radiance over
2014-2016, entirely before the 2017 production correction and long before the
2021 event. The same five brightness quintiles are then used for every annual
transition from 2015 through 2024.

No raster is re-extracted and no upstream production file is modified.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from diagnose_fixed_nuts2016_2021_reversal import load_locked_panel
from project_config import configured_path, load_config


BASELINE_YEARS = [2014, 2015, 2016]
Q_LABELS = [
    "Q1 darkest",
    "Q2",
    "Q3",
    "Q4",
    "Q5 brightest",
]
METRICS = {
    "raw": "radiance_mean_raw",
    "corrected": "radiance_mean_corrected",
}


def attach_baseline_brightness(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach one pre-2017 raw-radiance brightness reference to every region."""
    base = (
        panel[panel["year"].isin(BASELINE_YEARS)]
        .groupby(["NUTS_ID", "CNTR_CODE"], as_index=False)
        .agg(
            baseline_raw_nw=("radiance_mean_raw", "median"),
            baseline_n_years=("year", "nunique"),
        )
    )
    if not (base["baseline_n_years"] == len(BASELINE_YEARS)).all():
        raise RuntimeError("Baseline brightness is not complete for 2014-2016.")

    vals = pd.to_numeric(base["baseline_raw_nw"], errors="coerce")
    if (~np.isfinite(vals) | (vals <= 0)).any():
        raise RuntimeError("Baseline raw radiance contains non-positive values.")
    base["baseline_raw_nw"] = vals

    # qcut is used only to create five equally populated descriptive strata.
    # The continuous country-adjusted regression below does not depend on bins.
    try:
        base["brightness_quintile"] = pd.qcut(
            base["baseline_raw_nw"],
            q=5,
            labels=Q_LABELS,
            duplicates="raise",
        )
    except ValueError as exc:
        raise RuntimeError(
            "Could not form five distinct baseline-brightness quintiles."
        ) from exc

    if base["brightness_quintile"].isna().any():
        raise RuntimeError("At least one region lacks a brightness quintile.")

    base["log10_baseline_raw"] = np.log10(base["baseline_raw_nw"])
    centre = float(base["log10_baseline_raw"].median())
    base["log10_baseline_centered"] = base["log10_baseline_raw"] - centre

    merged = panel.merge(
        base[
            [
                "NUTS_ID",
                "baseline_raw_nw",
                "brightness_quintile",
                "log10_baseline_raw",
                "log10_baseline_centered",
            ]
        ],
        on="NUTS_ID",
        how="left",
        validate="many_to_one",
    )
    if merged["baseline_raw_nw"].isna().any():
        raise RuntimeError("Baseline brightness merge lost at least one region.")
    return merged, base


def make_transition_table(panel: pd.DataFrame) -> pd.DataFrame:
    d = panel.sort_values(["NUTS_ID", "year"]).copy()
    for label, col in METRICS.items():
        d[col] = pd.to_numeric(d[col], errors="coerce")
        if (~np.isfinite(d[col]) | (d[col] <= 0)).any():
            raise RuntimeError(f"Invalid values in {col}.")
        prev = d.groupby("NUTS_ID")[col].shift(1)
        d[f"previous_{label}_nw"] = prev
        d[f"delta_abs_{label}_nw"] = d[col] - prev

    d["previous_year"] = d.groupby("NUTS_ID")["year"].shift(1)
    d = d[d["previous_year"].notna()].copy()
    if ((d["year"].astype(int) - d["previous_year"].astype(int)) != 1).any():
        raise RuntimeError("Non-consecutive year pair in fixed panel.")

    expected_years = list(range(2015, 2025))
    observed = sorted(d["year"].astype(int).unique().tolist())
    if observed != expected_years:
        raise RuntimeError(
            f"Expected transition end years {expected_years}; found {observed}."
        )
    return d


def country_median_table(trans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, country), g in trans.groupby(["year", "CNTR_CODE"], observed=True):
        row = {
            "end_year": int(year),
            "CNTR_CODE": str(country),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
        }
        for metric in METRICS:
            x = g[f"delta_abs_{metric}_nw"]
            row[f"median_delta_{metric}_nw"] = float(x.median())
            row[f"mean_delta_{metric}_nw"] = float(x.mean())
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["end_year", "CNTR_CODE"])


def brightness_summary(trans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, quintile), g in trans.groupby(
        ["year", "brightness_quintile"], observed=True
    ):
        row_base = {
            "end_year": int(year),
            "brightness_quintile": str(quintile),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
            "baseline_min_nw": float(g["baseline_raw_nw"].min()),
            "baseline_median_nw": float(g["baseline_raw_nw"].median()),
            "baseline_max_nw": float(g["baseline_raw_nw"].max()),
        }
        for metric in METRICS:
            x = g[f"delta_abs_{metric}_nw"]
            row_base[f"median_delta_{metric}_nw"] = float(x.median())
            row_base[f"mean_delta_{metric}_nw"] = float(x.mean())
        rows.append(row_base)
    out = pd.DataFrame(rows)
    return out.sort_values(["end_year", "baseline_median_nw"])


def annual_summary(
    trans: pd.DataFrame,
    country: pd.DataFrame,
    by_brightness: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for year, g in trans.groupby("year", observed=True):
        row = {
            "end_year": int(year),
            "n_nuts3": int(g["NUTS_ID"].nunique()),
            "n_countries": int(g["CNTR_CODE"].nunique()),
        }
        c = country[country["end_year"] == int(year)]
        b = by_brightness[by_brightness["end_year"] == int(year)]

        for metric in METRICS:
            x = g[f"delta_abs_{metric}_nw"]
            cm = c[f"median_delta_{metric}_nw"]
            bm = b[f"median_delta_{metric}_nw"]
            row[f"median_delta_{metric}_nw"] = float(x.median())
            row[f"mean_delta_{metric}_nw"] = float(x.mean())
            row[f"pct_regions_increasing_{metric}"] = float(100.0 * (x > 0).mean())
            row[f"country_balanced_mean_median_{metric}_nw"] = float(cm.mean())
            row[f"country_median_of_medians_{metric}_nw"] = float(cm.median())
            row[f"brightness_quintile_median_range_{metric}_nw"] = float(
                bm.max() - bm.min()
            )

        rows.append(row)

    out = pd.DataFrame(rows).sort_values("end_year").reset_index(drop=True)
    for metric in METRICS:
        key = f"country_balanced_mean_median_{metric}_nw"
        out[f"abs_offset_rank_{metric}"] = (
            out[key].abs().rank(method="min", ascending=False).astype(int)
        )
        denom = out[key].abs()
        rng = out[f"brightness_quintile_median_range_{metric}_nw"]
        out[f"brightness_range_to_offset_ratio_{metric}"] = np.where(
            denom > 1e-9, rng / denom, np.nan
        )
    return out


def brightness_slope_table(trans: pd.DataFrame) -> pd.DataFrame:
    """Country-adjusted continuous brightness dependence for each transition."""
    rows = []
    for year, g in trans.groupby("year", observed=True):
        for metric in METRICS:
            response = f"delta_abs_{metric}_nw"
            fit = smf.ols(
                f"{response} ~ log10_baseline_centered + C(CNTR_CODE)",
                data=g,
            ).fit(cov_type="HC3")
            term = "log10_baseline_centered"
            ci = fit.conf_int().loc[term]
            rows.append(
                {
                    "end_year": int(year),
                    "metric": metric,
                    "n_nuts3": int(len(g)),
                    "slope_nw_per_log10_baseline": float(fit.params[term]),
                    "slope_se": float(fit.bse[term]),
                    "ci_low": float(ci.iloc[0]),
                    "ci_high": float(ci.iloc[1]),
                    "p": float(fit.pvalues[term]),
                    "r_squared": float(fit.rsquared),
                }
            )
    return pd.DataFrame(rows).sort_values(["end_year", "metric"])


def synthetic_self_test() -> None:
    """Exercise all downstream operations before reading project data."""
    rows = []
    # 100 regions guarantee five well-populated brightness quintiles.
    ids = [f"R{i:03d}" for i in range(100)]
    countries = np.repeat(["AA", "BB", "CC", "DD"], 25)
    baselines = np.geomspace(0.1, 20.0, 100)

    for nuts_id, country, baseline in zip(ids, countries, baselines):
        for year in range(2014, 2025):
            # Small proportional physical trend.
            physical = baseline * (1.0 + 0.005 * (year - 2014))

            # Synthetic known raw additive discontinuity from 2017 onward.
            raw = physical + (0.20 if year >= 2017 else 0.0)

            # Synthetic isolated 2021 additive pulse.
            if year == 2021:
                raw += 0.05

            # Corrected version removes only the known 2017 persistent offset.
            corrected = raw - (0.20 if year >= 2017 else 0.0)

            rows.append(
                {
                    "NUTS_ID": nuts_id,
                    "CNTR_CODE": country,
                    "year": year,
                    "radiance_mean_raw": raw,
                    "radiance_mean_corrected": corrected,
                }
            )

    panel = pd.DataFrame(rows)
    panel, _ = attach_baseline_brightness(panel)
    trans = make_transition_table(panel)
    country = country_median_table(trans)
    by_brightness = brightness_summary(trans)
    annual = annual_summary(trans, country, by_brightness)
    slopes = brightness_slope_table(trans)

    if len(trans) != 100 * 10:
        raise RuntimeError("Self-test failed: transition row count.")
    if len(by_brightness) != 10 * 5:
        raise RuntimeError("Self-test failed: brightness table row count.")
    if len(slopes) != 10 * 2:
        raise RuntimeError("Self-test failed: slope table row count.")

    a17 = annual.loc[annual["end_year"] == 2017].iloc[0]
    if not (
        abs(a17["country_balanced_mean_median_raw_nw"])
        > abs(a17["country_balanced_mean_median_corrected_nw"]) + 0.10
    ):
        raise RuntimeError("Self-test failed: 2017 positive-control attenuation.")

    a21 = annual.loc[annual["end_year"] == 2021].iloc[0]
    if abs(
        a21["country_balanced_mean_median_raw_nw"]
        - a21["country_balanced_mean_median_corrected_nw"]
    ) > 1e-10:
        raise RuntimeError("Self-test failed: persistent correction should cancel in 2021.")

    print("Absolute-offset downstream synthetic self-test: PASS", flush=True)


def make_plots(
    trans: pd.DataFrame,
    annual: pd.DataFrame,
    by_brightness: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    ax.plot(
        annual["end_year"],
        annual["country_balanced_mean_median_raw_nw"],
        marker="o",
        label="raw average_masked",
    )
    ax.plot(
        annual["end_year"],
        annual["country_balanced_mean_median_corrected_nw"],
        marker="o",
        label="production-corrected",
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(2017, linestyle="--", linewidth=0.9)
    ax.axvline(2021, linestyle="--", linewidth=0.9)
    ax.set_xlabel("Transition ending in year")
    ax.set_ylabel("Country-balanced median absolute change (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Annual absolute VIIRS change on fixed NUTS2016 geography")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "annual_absolute_change_raw_vs_corrected.pdf", dpi=dpi)
    plt.close(fig)

    for metric in ["raw", "corrected"]:
        fig, ax = plt.subplots(figsize=(9.0, 6.0))
        for label in Q_LABELS:
            g = by_brightness[
                by_brightness["brightness_quintile"].astype(str) == label
            ]
            ax.plot(
                g["end_year"],
                g[f"median_delta_{metric}_nw"],
                marker="o",
                label=label,
            )
        ax.axhline(0, linewidth=0.8)
        ax.axvline(2017, linestyle="--", linewidth=0.9)
        ax.axvline(2021, linestyle="--", linewidth=0.9)
        ax.set_xlabel("Transition ending in year")
        ax.set_ylabel("Median absolute change (nW cm$^{-2}$ sr$^{-1}$)")
        ax.set_title(
            f"{metric.capitalize()} absolute change by fixed 2014-2016 brightness quintile"
        )
        ax.legend()
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(
            output_dir / f"{metric}_absolute_change_by_baseline_brightness.pdf",
            dpi=dpi,
        )
        plt.close(fig)

    g = trans[trans["year"] == 2021].copy()
    fig, ax = plt.subplots(figsize=(7.6, 6.3))
    ax.scatter(
        g["baseline_raw_nw"],
        g["delta_abs_raw_nw"],
        s=12,
        alpha=0.4,
    )
    med = (
        g.groupby("brightness_quintile", observed=True)
        .agg(
            baseline=("baseline_raw_nw", "median"),
            delta=("delta_abs_raw_nw", "median"),
        )
        .reset_index()
        .sort_values("baseline")
    )
    ax.plot(med["baseline"], med["delta"], marker="o", linewidth=1.5)
    ax.axhline(0, linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Fixed 2014-2016 median raw radiance (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_ylabel("Raw absolute change 2020→2021 (nW cm$^{-2}$ sr$^{-1}$)")
    ax.set_title("Is the 2021 raw-radiance shift approximately additive?")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "raw_2021_absolute_change_vs_baseline_brightness.pdf", dpi=dpi)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    # Fail cheaply if downstream pandas/statsmodels transformations are broken.
    synthetic_self_test()

    config = load_config(args.config)
    diag_cfg = config["analysis"]["fixed_nuts2016_2021_diagnostic"]
    panel = load_locked_panel(config)

    panel, baseline = attach_baseline_brightness(panel)
    trans = make_transition_table(panel)
    country = country_median_table(trans)
    by_brightness = brightness_summary(trans)
    annual = annual_summary(trans, country, by_brightness)
    slopes = brightness_slope_table(trans)

    expected_n = int(diag_cfg["expected_n_nuts3"])
    if trans["NUTS_ID"].nunique() != expected_n:
        raise RuntimeError(
            f"Expected {expected_n} NUTS3 regions; found {trans['NUTS_ID'].nunique()}."
        )
    if len(trans) != expected_n * 10:
        raise RuntimeError(
            f"Expected {expected_n * 10} transition rows; found {len(trans)}."
        )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_absolute_offset_audit"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    outputs = {
        "baseline": processed / "fixed_nuts2016_baseline_brightness_2014_2016.csv",
        "transitions": processed
        / "fixed_nuts2016_absolute_transition_changes_2015_2024.csv",
        "annual": processed
        / "fixed_nuts2016_absolute_change_annual_summary_2015_2024.csv",
        "brightness": processed
        / "fixed_nuts2016_absolute_change_by_brightness_2015_2024.csv",
        "country": processed
        / "fixed_nuts2016_absolute_change_country_medians_2015_2024.csv",
        "slopes": processed
        / "fixed_nuts2016_absolute_change_brightness_slopes_2015_2024.csv",
        "summary": processed / "fixed_nuts2016_absolute_offset_audit_summary.txt",
    }

    baseline.to_csv(outputs["baseline"], index=False)
    trans.to_csv(outputs["transitions"], index=False)
    annual.to_csv(outputs["annual"], index=False)
    by_brightness.to_csv(outputs["brightness"], index=False)
    country.to_csv(outputs["country"], index=False)
    slopes.to_csv(outputs["slopes"], index=False)

    make_plots(
        trans,
        annual,
        by_brightness,
        figures,
        int(diag_cfg["figure_dpi"]),
    )

    a17 = annual.loc[annual["end_year"] == 2017].iloc[0]
    a21 = annual.loc[annual["end_year"] == 2021].iloc[0]
    s21_raw = slopes[
        (slopes["end_year"] == 2021) & (slopes["metric"] == "raw")
    ].iloc[0]
    b21 = by_brightness[by_brightness["end_year"] == 2021].copy()

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 ABSOLUTE-OFFSET AUDIT: 2015-2024\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Locked panel: {trans['NUTS_ID'].nunique()} NUTS3 regions; "
            f"{trans['CNTR_CODE'].nunique()} countries.\n"
        )
        h.write(
            "Primary response: unmodified raw EOG average_masked radiance.\n"
        )
        h.write(
            "Fixed brightness reference: median raw radiance over 2014-2016.\n"
        )
        h.write(
            "Purpose: test whether 2020->2021 is an unusually large, "
            "approximately additive absolute radiance shift across baseline "
            "brightness.\n\n"
        )

        h.write("ANNUAL ABSOLUTE-CHANGE SUMMARY\n")
        h.write(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nCOUNTRY-ADJUSTED BRIGHTNESS SLOPES\n")
        h.write(slopes.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\n2021 BY BASELINE-BRIGHTNESS QUINTILE\n")
        h.write(
            b21.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )

        h.write("\n\nPRE-SPECIFIED CHECKS\n")
        h.write(
            f"- 2017 raw country-balanced median shift: "
            f"{a17['country_balanced_mean_median_raw_nw']:+.6f} nW.\n"
        )
        h.write(
            f"- 2017 corrected country-balanced median shift: "
            f"{a17['country_balanced_mean_median_corrected_nw']:+.6f} nW.\n"
        )
        h.write(
            f"- 2021 raw country-balanced median shift: "
            f"{a21['country_balanced_mean_median_raw_nw']:+.6f} nW.\n"
        )
        h.write(
            f"- 2021 raw absolute-offset rank among 2015-2024 transitions: "
            f"{int(a21['abs_offset_rank_raw'])} (1 = largest absolute shift).\n"
        )
        h.write(
            f"- 2021 range of raw quintile median shifts: "
            f"{a21['brightness_quintile_median_range_raw_nw']:.6f} nW.\n"
        )
        h.write(
            f"- 2021 country-adjusted raw brightness slope: "
            f"{s21_raw['slope_nw_per_log10_baseline']:+.6f} nW per log10 unit "
            f"(95% CI {s21_raw['ci_low']:+.6f} to {s21_raw['ci_high']:+.6f}).\n"
        )

        h.write(
            "\nDecision rule:\n"
            "- The known 2017 raw discontinuity is a positive control: the "
            "production-corrected series should attenuate it.\n"
            "- A 2021 shift that is large in absolute magnitude but has a "
            "small brightness slope and a small range across fixed brightness "
            "quintiles supports an additive annual-product disturbance.\n"
            "- If comparable near-flat absolute shifts occur in several other "
            "years, treat them as normal annual-product variability rather than "
            "correcting 2021 alone.\n"
            "- Only after this audit should we decide whether a specific "
            "product correction, exclusion/sensitivity analysis, or a model "
            "that explicitly absorbs annual measurement shocks is justified.\n"
        )

    print("\nFIXED NUTS2016 ABSOLUTE-OFFSET AUDIT")
    print(
        f"Locked panel: {trans['NUTS_ID'].nunique()} NUTS3; "
        f"{trans['CNTR_CODE'].nunique()} countries"
    )
    print("\nANNUAL ABSOLUTE-CHANGE SUMMARY")
    print(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nCOUNTRY-ADJUSTED BRIGHTNESS SLOPES")
    print(slopes.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\n2021 BY BASELINE-BRIGHTNESS QUINTILE")
    print(b21.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
