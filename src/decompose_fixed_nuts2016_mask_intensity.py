#!/usr/bin/env python3
"""Decompose annual raw VIIRS mean change into mask occupancy and lit-pixel intensity.

This stage uses only existing fixed-NUTS2016 outputs. It does NOT re-extract
rasters and does NOT modify production data.

For each NUTS3-year, the raw average_masked mean R can be written exactly as

    R = p * I

where p is the physical-area-weighted fraction of the polygon with
average_masked radiance > 0 and I is the corresponding conditional mean
radiance over positive area.

For adjacent years 0 and 1, the exact midpoint/Shapley decomposition is

    delta_R_occupancy = (p1 - p0) * (I0 + I1) / 2
    delta_R_intensity = (I1 - I0) * (p0 + p1) / 2

and

    delta_R = delta_R_occupancy + delta_R_intensity.

This provides a direct, scale-preserving test of how much of the 2021 pulse is
associated with changing masked-positive area versus radiance changes within
already-positive area.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from diagnose_fixed_nuts2016_2021_reversal import load_locked_panel
from diagnose_fixed_nuts2016_mask_occupancy import attach_fixed_brightness
from project_config import configured_path, load_config


OCCUPANCY_FILE = "fixed_nuts2016_mask_positive_area_2014_2024.csv"


def load_inputs(config: dict) -> pd.DataFrame:
    processed = configured_path(config, "processed")
    occ_path = processed / OCCUPANCY_FILE
    if not occ_path.exists():
        raise FileNotFoundError(
            f"Missing {occ_path}. Run "
            "'python src/run_pipeline.py future-radiance-mask-occupancy-audit' first."
        )

    panel = attach_fixed_brightness(load_locked_panel(config))
    occ = pd.read_csv(occ_path)
    required = {"year", "NUTS_ID", "positive_masked_area_fraction"}
    missing = required - set(occ.columns)
    if missing:
        raise RuntimeError(
            f"Mask occupancy file is missing columns: {sorted(missing)}"
        )

    occ = occ[list(required)].copy()
    occ["NUTS_ID"] = occ["NUTS_ID"].astype(str)
    panel["NUTS_ID"] = panel["NUTS_ID"].astype(str)

    if occ.duplicated(["year", "NUTS_ID"]).any():
        raise RuntimeError("Duplicate year-NUTS_ID rows in mask occupancy file.")

    d = panel.merge(
        occ,
        on=["year", "NUTS_ID"],
        how="left",
        validate="one_to_one",
    )
    if d["positive_masked_area_fraction"].isna().any():
        raise RuntimeError("Mask occupancy merge lost at least one panel row.")

    p = pd.to_numeric(d["positive_masked_area_fraction"], errors="coerce")
    r = pd.to_numeric(d["radiance_mean_raw"], errors="coerce")
    if (~np.isfinite(p) | (p <= 0) | (p > 1 + 1e-8)).any():
        raise RuntimeError("Invalid positive masked-area fraction.")
    if (~np.isfinite(r) | (r <= 0)).any():
        raise RuntimeError("Invalid raw zonal radiance.")

    d["positive_masked_area_fraction"] = p
    d["radiance_mean_raw"] = r
    d["positive_area_mean_raw_nw"] = r / p

    # Algebraic reconstruction is an invariant. This is intentionally strict:
    # the decomposition is meaningful only if R = p * I on the same support.
    reconstructed = p * d["positive_area_mean_raw_nw"]
    err = np.abs(reconstructed - r)
    if float(err.max()) > 1e-10:
        raise RuntimeError(
            f"R = p*I reconstruction failed; max absolute error={err.max():.3g}"
        )

    return d.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)


def make_transitions(d: pd.DataFrame) -> pd.DataFrame:
    x = d.copy()
    for col in [
        "radiance_mean_raw",
        "positive_masked_area_fraction",
        "positive_area_mean_raw_nw",
    ]:
        x[f"previous_{col}"] = x.groupby("NUTS_ID")[col].shift(1)

    x["previous_year"] = x.groupby("NUTS_ID")["year"].shift(1)
    t = x[x["previous_year"].notna()].copy()
    if ((t["year"].astype(int) - t["previous_year"].astype(int)) != 1).any():
        raise RuntimeError("Non-consecutive annual transition.")

    p0 = t["previous_positive_masked_area_fraction"]
    p1 = t["positive_masked_area_fraction"]
    i0 = t["previous_positive_area_mean_raw_nw"]
    i1 = t["positive_area_mean_raw_nw"]

    t["delta_raw_nw"] = (
        t["radiance_mean_raw"] - t["previous_radiance_mean_raw"]
    )
    t["occupancy_component_nw"] = (p1 - p0) * (i0 + i1) / 2.0
    t["intensity_component_nw"] = (i1 - i0) * (p0 + p1) / 2.0
    t["decomposition_residual_nw"] = (
        t["delta_raw_nw"]
        - t["occupancy_component_nw"]
        - t["intensity_component_nw"]
    )

    max_resid = float(np.abs(t["decomposition_residual_nw"]).max())
    if max_resid > 1e-9:
        raise RuntimeError(
            f"Exact decomposition invariant failed; max residual={max_resid:.3g}"
        )

    return t


def country_table(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, country), g in t.groupby(["year", "CNTR_CODE"], observed=True):
        rows.append(
            {
                "end_year": int(year),
                "CNTR_CODE": str(country),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "median_total_nw": float(g["delta_raw_nw"].median()),
                "median_occupancy_nw": float(
                    g["occupancy_component_nw"].median()
                ),
                "median_intensity_nw": float(
                    g["intensity_component_nw"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def annual_summary(t: pd.DataFrame, countries: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for year, g in t.groupby("year", observed=True):
        c = countries[countries["end_year"] == int(year)]
        occ_cb = float(c["median_occupancy_nw"].mean())
        int_cb = float(c["median_intensity_nw"].mean())
        total_cb = occ_cb + int_cb
        rows.append(
            {
                "end_year": int(year),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "median_total_nw": float(g["delta_raw_nw"].median()),
                "median_occupancy_component_nw": float(
                    g["occupancy_component_nw"].median()
                ),
                "median_intensity_component_nw": float(
                    g["intensity_component_nw"].median()
                ),
                "country_balanced_occupancy_component_nw": occ_cb,
                "country_balanced_intensity_component_nw": int_cb,
                "country_balanced_total_from_components_nw": total_cb,
                "occupancy_fraction_of_cb_total": (
                    occ_cb / total_cb if abs(total_cb) > 1e-12 else np.nan
                ),
                "intensity_fraction_of_cb_total": (
                    int_cb / total_cb if abs(total_cb) > 1e-12 else np.nan
                ),
                "median_abs_decomposition_residual_nw": float(
                    np.abs(g["decomposition_residual_nw"]).median()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("end_year").reset_index(drop=True)


def grouped_summary(t: pd.DataFrame, group: str) -> pd.DataFrame:
    rows = []
    for (year, label), g in t.groupby(["year", group], observed=True):
        occ = float(g["occupancy_component_nw"].median())
        inten = float(g["intensity_component_nw"].median())
        total_components = occ + inten
        rows.append(
            {
                "end_year": int(year),
                group: str(label),
                "n_nuts3": int(g["NUTS_ID"].nunique()),
                "median_total_nw": float(g["delta_raw_nw"].median()),
                "median_occupancy_component_nw": occ,
                "median_intensity_component_nw": inten,
                "occupancy_fraction_of_median_components": (
                    occ / total_components
                    if abs(total_components) > 1e-12
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def synthetic_self_test() -> None:
    rows = []
    # Construct exact R=p*I data with one occupancy-only and one intensity-only
    # event so both decomposition terms are checked independently.
    for i in range(20):
        nuts = f"R{i:02d}"
        country = "AA" if i < 10 else "BB"
        brightness = 0.2 + 0.1 * i
        quint = ["Q1 darkest", "Q2", "Q3", "Q4", "Q5 brightest"][min(i // 4, 4)]
        density = "<30" if i < 10 else "30-100"
        for year in range(2014, 2025):
            p = 0.20 + i * 0.002
            intensity = brightness
            if year == 2021:
                p += 0.02
            if year == 2022:
                p += 0.02
                intensity += 0.10
            rows.append(
                {
                    "NUTS_ID": nuts,
                    "CNTR_CODE": country,
                    "year": year,
                    "density_band": density,
                    "brightness_quintile": quint,
                    "baseline_raw_nw": brightness,
                    "positive_masked_area_fraction": p,
                    "positive_area_mean_raw_nw": intensity,
                    "radiance_mean_raw": p * intensity,
                }
            )

    d = pd.DataFrame(rows)
    t = make_transitions(d)

    g21 = t[t["year"] == 2021]
    if not (g21["occupancy_component_nw"] > 0).all():
        raise RuntimeError("Self-test failed to recover occupancy-only pulse.")
    if float(np.abs(g21["intensity_component_nw"]).max()) > 1e-12:
        raise RuntimeError("Self-test produced false intensity change in 2021.")

    g22 = t[t["year"] == 2022]
    if not (g22["intensity_component_nw"] > 0).all():
        raise RuntimeError("Self-test failed to recover intensity-only pulse.")
    if float(np.abs(g22["occupancy_component_nw"]).max()) > 1e-12:
        raise RuntimeError("Self-test produced false occupancy change in 2022.")

    if float(np.abs(t["decomposition_residual_nw"]).max()) > 1e-12:
        raise RuntimeError("Self-test decomposition does not sum exactly.")

    countries = country_table(t)
    annual = annual_summary(t, countries)
    by_brightness = grouped_summary(t, "brightness_quintile")
    by_density = grouped_summary(t, "density_band")
    if annual.empty or by_brightness.empty or by_density.empty:
        raise RuntimeError("Self-test summary generation failed.")

    print("Mask/intensity decomposition synthetic self-test: PASS", flush=True)


def make_plots(
    annual: pd.DataFrame,
    by_brightness: pd.DataFrame,
    by_density: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.2, 5.9))
    ax.plot(
        annual["end_year"],
        annual["country_balanced_occupancy_component_nw"],
        marker="o",
        label="mask occupancy component",
    )
    ax.plot(
        annual["end_year"],
        annual["country_balanced_intensity_component_nw"],
        marker="o",
        label="within-positive intensity component",
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(2021, linestyle="--", linewidth=0.9)
    ax.set_xlabel("Transition ending in year")
    ax.set_ylabel(
        "Country-balanced absolute-change component "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title("Exact decomposition of annual raw VIIRS mean change")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "annual_occupancy_vs_intensity_components.pdf", dpi=dpi)
    plt.close(fig)

    for table, group, filename, xlabel in [
        (
            by_brightness,
            "brightness_quintile",
            "decomposition_2021_by_brightness.pdf",
            "Fixed 2014-2016 brightness quintile",
        ),
        (
            by_density,
            "density_band",
            "decomposition_2021_by_density.pdf",
            "Fixed population-density band (people km$^{-2}$)",
        ),
    ]:
        g = table[table["end_year"] == 2021].copy()
        if group == "brightness_quintile":
            order = ["Q1 darkest", "Q2", "Q3", "Q4", "Q5 brightest"]
        else:
            order = ["<30", "30-100", "100-300", "300-1000", ">=1000"]
        g[group] = pd.Categorical(g[group], categories=order, ordered=True)
        g = g.sort_values(group)

        x = np.arange(len(g))
        width = 0.36
        fig, ax = plt.subplots(figsize=(9.0, 5.9))
        ax.bar(
            x - width / 2,
            g["median_occupancy_component_nw"],
            width=width,
            label="mask occupancy",
        )
        ax.bar(
            x + width / 2,
            g["median_intensity_component_nw"],
            width=width,
            label="within-positive intensity",
        )
        ax.axhline(0, linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(g[group].astype(str))
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Median 2020→2021 component (nW cm$^{-2}$ sr$^{-1}$)")
        ax.set_title("2021 raw-radiance change decomposition")
        ax.legend()
        ax.grid(alpha=0.2, axis="y")
        fig.tight_layout()
        fig.savefig(output_dir / filename, dpi=dpi)
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    synthetic_self_test()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_2021_diagnostic"]

    d = load_inputs(config)
    t = make_transitions(d)
    countries = country_table(t)
    annual = annual_summary(t, countries)
    by_brightness = grouped_summary(t, "brightness_quintile")
    by_density = grouped_summary(t, "density_band")

    expected_n = int(cfg["expected_n_nuts3"])
    if t["NUTS_ID"].nunique() != expected_n:
        raise RuntimeError(
            f"Expected {expected_n} fixed NUTS3 regions; "
            f"found {t['NUTS_ID'].nunique()}."
        )
    if len(t) != expected_n * 10:
        raise RuntimeError(
            f"Expected {expected_n * 10} transition rows; found {len(t)}."
        )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_mask_intensity_decomposition"
    )
    processed.mkdir(parents=True, exist_ok=True)

    outputs = {
        "regional": processed
        / "fixed_nuts2016_mask_intensity_decomposition_2015_2024.csv",
        "annual": processed
        / "fixed_nuts2016_mask_intensity_decomposition_annual_2015_2024.csv",
        "country": processed
        / "fixed_nuts2016_mask_intensity_decomposition_country_2015_2024.csv",
        "brightness": processed
        / "fixed_nuts2016_mask_intensity_decomposition_by_brightness_2015_2024.csv",
        "density": processed
        / "fixed_nuts2016_mask_intensity_decomposition_by_density_2015_2024.csv",
        "summary": processed
        / "fixed_nuts2016_mask_intensity_decomposition_summary.txt",
    }

    t.to_csv(outputs["regional"], index=False)
    annual.to_csv(outputs["annual"], index=False)
    countries.to_csv(outputs["country"], index=False)
    by_brightness.to_csv(outputs["brightness"], index=False)
    by_density.to_csv(outputs["density"], index=False)

    make_plots(
        annual,
        by_brightness,
        by_density,
        figures,
        int(cfg["figure_dpi"]),
    )

    a21 = annual.loc[annual["end_year"] == 2021].iloc[0]
    a22 = annual.loc[annual["end_year"] == 2022].iloc[0]
    b21 = by_brightness[by_brightness["end_year"] == 2021]
    d21 = by_density[by_density["end_year"] == 2021]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 MASK/INTENSITY DECOMPOSITION: 2015-2024\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Locked panel: {t['NUTS_ID'].nunique()} NUTS3 regions; "
            f"{t['CNTR_CODE'].nunique()} countries.\n"
        )
        h.write(
            "Identity: raw average_masked zonal mean = positive masked-area "
            "fraction × conditional mean radiance over positive area.\n"
        )
        h.write(
            "Annual change is partitioned exactly with a symmetric midpoint/"
            "Shapley decomposition into mask-occupancy and within-positive "
            "intensity components.\n\n"
        )

        h.write("ANNUAL COMPONENT SUMMARY\n")
        h.write(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\n2021 BY FIXED BASELINE BRIGHTNESS\n")
        h.write(b21.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\n2021 BY FIXED POPULATION-DENSITY BAND\n")
        h.write(d21.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

        h.write("\n\nPRE-SPECIFIED 2021/2022 CHECKS\n")
        h.write(
            f"- 2021 country-balanced occupancy component: "
            f"{a21['country_balanced_occupancy_component_nw']:+.6f} nW.\n"
        )
        h.write(
            f"- 2021 country-balanced intensity component: "
            f"{a21['country_balanced_intensity_component_nw']:+.6f} nW.\n"
        )
        h.write(
            f"- 2021 occupancy fraction of country-balanced decomposed total: "
            f"{100*a21['occupancy_fraction_of_cb_total']:+.2f}%.\n"
        )
        h.write(
            f"- 2022 country-balanced occupancy component: "
            f"{a22['country_balanced_occupancy_component_nw']:+.6f} nW.\n"
        )
        h.write(
            f"- 2022 country-balanced intensity component: "
            f"{a22['country_balanced_intensity_component_nw']:+.6f} nW.\n"
        )
        h.write(
            f"- 2022 occupancy fraction of country-balanced decomposed total: "
            f"{100*a22['occupancy_fraction_of_cb_total']:+.2f}%.\n"
        )

        h.write(
            "\nDecision rule:\n"
            "- A dominant 2021 occupancy component would show that the annual "
            "background/masked-positive support itself generates most of the "
            "pulse.\n"
            "- A dominant intensity component would show that most of the "
            "pulse occurs within area already retained as positive.\n"
            "- Mixed components imply that both mask support and retained-pixel "
            "radiance changed; any correction must then target the actual "
            "processing mechanism rather than subtracting one arbitrary NUTS3 "
            "offset.\n"
        )

    print("\nFIXED NUTS2016 MASK/INTENSITY DECOMPOSITION")
    print(
        f"Locked panel: {t['NUTS_ID'].nunique()} NUTS3; "
        f"{t['CNTR_CODE'].nunique()} countries"
    )
    print("\nANNUAL COMPONENT SUMMARY")
    print(annual.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\n2021 BY FIXED BASELINE BRIGHTNESS")
    print(b21.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\n2021 BY FIXED POPULATION-DENSITY BAND")
    print(d21.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
