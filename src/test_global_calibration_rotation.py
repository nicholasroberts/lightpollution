#!/usr/bin/env python3
"""Test whether the density-dependent VIIRS rotation survives global 2017 correction.

IMPORTANT
---------
This is a diagnostic stage. It does NOT overwrite the production calibrated
radiance table.

The EOG annual VIIRS series has a known 2017 calibration discontinuity that is
not country-specific. The current production branch historically corrected only
five comparison countries; this stage explicitly tests the consequences.

Five radiance scenarios are compared on identical population data:

1. raw_uncorrected
2. legacy_five_country (the current historical production column)
3. published_global_0p15: subtract 0.15 nW cm-2 sr-1 from every NUTS3 mean
   from 2017 onward
4. empirical_global_dark20: one Europe-wide offset estimated from the darkest
   20% of matched 2015-2018 NUTS3 regions
5. empirical_country_median_dark20: estimate the same dark-region offset within
   every country with adequate data, take the median country offset, then apply
   that ONE common offset to all regions from 2017 onward

The global corrections are intentionally applied at NUTS3-mean level in this
diagnostic because that matches the existing pipeline and directly answers
whether the previously reported population-density rotation is robust to a
global additive correction. The appropriate production correction for the
masked annual raster will be chosen only after inspecting these results.

For every scenario the script reports:
- growth (% yr-1) in the established population-density bands;
- the band x year likelihood-ratio test;
- a continuous cubic density x year likelihood-ratio test;
- modelled growth at reference population densities;
- annual country-adjusted log-log population-radiance slope.

The analysis is run both for the full 2013-2024 record and for fixed NUTS2016
geography over 2013-2020.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


CAL_YEARS = [2015, 2016, 2017, 2018]


def drop_unused_categories(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    for col in out.columns:
        if isinstance(out[col].dtype, pd.CategoricalDtype):
            out[col] = out[col].cat.remove_unused_categories()
    return out


def fit_mixedlm(formula: str, data: pd.DataFrame):
    fit = drop_unused_categories(data)
    model = smf.mixedlm(
        formula,
        data=fit,
        groups=fit["NUTS_ID"],
        re_formula="1",
    )
    last = None
    for method in ("lbfgs", "bfgs", "cg"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = model.fit(
                    reml=False,
                    method=method,
                    maxiter=3500,
                    disp=False,
                )
            if result.converged:
                return result, method
            last = RuntimeError(f"non-converged with {method}")
        except Exception as exc:
            last = exc
    raise RuntimeError(f"MixedLM failed for {formula}: {last}")


def lrt(reduced, full):
    lr = max(0.0, 2.0 * (full.llf - reduced.llf))
    df = int(len(full.fe_params) - len(reduced.fe_params))
    p = float(chi2.sf(lr, df)) if df > 0 else np.nan
    return lr, df, p


def contrast(result, weights):
    names = list(result.fe_params.index)
    L = np.zeros(len(names), dtype=float)
    for name, value in weights.items():
        if name not in names:
            raise KeyError(f"Missing parameter {name!r}")
        L[names.index(name)] = value
    beta = result.fe_params.to_numpy(float)
    cov = result.cov_params().loc[names, names].to_numpy(float)
    est = float(L @ beta)
    var = float(L @ cov @ L)
    se = float(np.sqrt(max(0.0, var)))
    z = est / se if se > 0 else np.nan
    p = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return est, se, p


def candidate_offset(wide: pd.DataFrame) -> float:
    d15 = (wide[2016] - wide[2015]).median()
    d16 = (wide[2017] - wide[2016]).median()
    d17 = (wide[2018] - wide[2017]).median()
    expected = np.nanmean([d15, d17])
    return float(d16 - expected)


def estimate_global_offsets(config: dict):
    cfg = config["analysis"]["global_calibration_rotation_test"]
    raw_path = (
        configured_path(config, "processed")
        / config["outputs"]["radiance_by_nuts3"]
    )
    raw = pd.read_parquet(raw_path)
    cal = raw[raw["year"].isin(CAL_YEARS)].copy()

    wide = (
        cal.pivot_table(
            index=["NUTS_ID", "CNTR_CODE"],
            columns="year",
            values="radiance_mean",
            aggfunc="first",
        )
        .dropna(subset=CAL_YEARS)
        .reset_index()
    )

    dark_fraction = float(cfg["dark_fraction"])
    threshold = float(wide[2016].quantile(dark_fraction))
    dark = wide[wide[2016] <= threshold].copy()
    global_offset = candidate_offset(dark)

    min_country = int(cfg["min_country_matched_regions"])
    country_rows = []
    for country, g in wide.groupby("CNTR_CODE"):
        if len(g) < min_country:
            continue
        t = float(g[2016].quantile(dark_fraction))
        gd = g[g[2016] <= t].copy()
        if len(gd) < 2:
            continue
        country_rows.append(
            {
                "CNTR_CODE": country,
                "n_matched": len(g),
                "n_dark": len(gd),
                "dark_threshold_2016": t,
                "candidate_offset": candidate_offset(gd),
            }
        )
    countries = pd.DataFrame(country_rows).sort_values("candidate_offset")
    median_country_offset = float(countries["candidate_offset"].median())

    summary = pd.DataFrame(
        [
            {
                "scenario": "published_global_0p15",
                "offset_nW_cm2_sr": float(cfg["published_global_offset"]),
                "estimation": "published literature sensitivity value",
                "n_matched": np.nan,
                "n_dark": np.nan,
            },
            {
                "scenario": "empirical_global_dark20",
                "offset_nW_cm2_sr": global_offset,
                "estimation": "darkest fraction of all matched European NUTS3",
                "n_matched": len(wide),
                "n_dark": len(dark),
            },
            {
                "scenario": "empirical_country_median_dark20",
                "offset_nW_cm2_sr": median_country_offset,
                "estimation": "median of within-country dark-region offsets",
                "n_matched": int(countries["n_matched"].sum()),
                "n_dark": int(countries["n_dark"].sum()),
            },
        ]
    )
    return summary, countries


def add_scenarios(merged: pd.DataFrame, offsets: pd.DataFrame, config: dict):
    d = merged.copy()
    d["rad_raw_uncorrected"] = d["radiance_mean_raw"].astype(float)
    d["rad_legacy_five_country"] = d["radiance_mean_corrected"].astype(float)

    off = dict(zip(offsets["scenario"], offsets["offset_nW_cm2_sr"]))
    post = d["year"] >= int(config["calibration"]["apply_from_year"])

    for scenario in [
        "published_global_0p15",
        "empirical_global_dark20",
        "empirical_country_median_dark20",
    ]:
        col = f"rad_{scenario}"
        d[col] = d["radiance_mean_raw"].astype(float)
        d.loc[post, col] = d.loc[post, "radiance_mean_raw"].astype(float) - off[scenario]

    return d


def prepare_scenario(
    data: pd.DataFrame,
    rad_col: str,
    config: dict,
    first_year: int,
    last_year: int,
    nuts_release: int | None,
    common_positive_mask: pd.Series | None = None,
):
    cfg = config["analysis"]["global_calibration_rotation_test"]
    d = data[data["year"].between(first_year, last_year, inclusive="both")].copy()
    if nuts_release is not None:
        d = d[d["nuts_release"] == nuts_release].copy()

    if common_positive_mask is not None:
        d = d.loc[common_positive_mask.reindex(d.index).fillna(False)].copy()

    ok = (
        np.isfinite(pd.to_numeric(d[rad_col], errors="coerce"))
        & np.isfinite(pd.to_numeric(d["population_density"], errors="coerce"))
        & (d[rad_col] > 0)
        & (d["population_density"] > 0)
        & d["NUTS_ID"].notna()
        & d["CNTR_CODE"].notna()
    )
    d = d.loc[ok].copy()

    repeated = d.groupby("NUTS_ID").size()
    d = d[
        d["NUTS_ID"].isin(
            repeated[repeated >= int(cfg["min_repeated_years"])].index
        )
    ].copy()
    counts = d.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    d = d[
        d["CNTR_CODE"].isin(
            counts[counts >= int(cfg["min_country_regions"])].index
        )
    ].copy()

    d["log_rad"] = np.log10(d[rad_col].astype(float))
    d["log_pop"] = np.log10(d["population_density"].astype(float))
    d["year_c"] = d["year"].astype(float) - ((first_year + last_year) / 2.0)
    d["country"] = d["CNTR_CODE"].astype("category")

    ref = (
        d.groupby("NUTS_ID", as_index=False)["population_density"]
        .median()
        .rename(columns={"population_density": "reference_population_density"})
    )
    d = d.merge(ref, on="NUTS_ID", how="left", validate="many_to_one")

    edges = [-np.inf] + list(cfg["density_band_breaks"]) + [np.inf]
    d["density_band"] = pd.cut(
        d["reference_population_density"],
        bins=edges,
        labels=cfg["density_band_labels"],
        right=False,
        ordered=True,
    )

    center = float(d["log_pop"].median())
    d["log_pop_c"] = d["log_pop"] - center
    d.attrs["log_pop_center"] = center
    return d


def fit_band_model(d: pd.DataFrame, labels):
    reduced, _ = fit_mixedlm(
        "log_rad ~ year_c + C(density_band) + C(country)",
        d,
    )
    full, _ = fit_mixedlm(
        "log_rad ~ year_c * C(density_band) + C(country)",
        d,
    )
    lr, df, p = lrt(reduced, full)

    rows = []
    ref = labels[0]
    for band in labels:
        weights = {"year_c": 1.0}
        if band != ref:
            candidates = [
                f"year_c:C(density_band)[T.{band}]",
                f"C(density_band)[T.{band}]:year_c",
            ]
            name = next(
                (x for x in candidates if x in full.fe_params.index),
                None,
            )
            if name is None:
                continue
            weights[name] = 1.0
        est, se, pp = contrast(full, weights)
        rows.append(
            {
                "density_band": band,
                "log10_radiance_trend_per_year": est,
                "trend_se": se,
                "trend_p": pp,
                "pct_radiance_change_per_year": 100.0 * (10.0**est - 1.0),
            }
        )

    return pd.DataFrame(rows), {
        "band_year_lr": lr,
        "band_year_df": df,
        "band_year_p": p,
    }


def fit_continuous_model(d: pd.DataFrame, reference_densities):
    base = (
        "log_rad ~ C(country) + year_c + log_pop_c "
        "+ I(log_pop_c ** 2) + I(log_pop_c ** 3)"
    )
    full = (
        "log_rad ~ C(country) + "
        "year_c * (log_pop_c + I(log_pop_c ** 2) + I(log_pop_c ** 3))"
    )
    m0, _ = fit_mixedlm(base, d)
    m1, _ = fit_mixedlm(full, d)
    lr, df, p = lrt(m0, m1)

    center = d.attrs["log_pop_center"]
    rows = []
    for density in reference_densities:
        x = np.log10(float(density)) - center
        weights = {
            "year_c": 1.0,
            "year_c:log_pop_c": x,
            "year_c:I(log_pop_c ** 2)": x**2,
            "year_c:I(log_pop_c ** 3)": x**3,
        }
        est, se, pp = contrast(m1, weights)
        rows.append(
            {
                "population_density": float(density),
                "log10_radiance_trend_per_year": est,
                "trend_se": se,
                "trend_p": pp,
                "pct_radiance_change_per_year": 100.0 * (10.0**est - 1.0),
            }
        )
    return pd.DataFrame(rows), {
        "continuous_density_year_lr": lr,
        "continuous_density_year_df": df,
        "continuous_density_year_p": p,
    }


def annual_common_slopes(d: pd.DataFrame):
    rows = []
    for year, g in d.groupby("year"):
        fit = drop_unused_categories(g)
        if fit["country"].nunique() < 2:
            continue
        result = smf.ols(
            "log_rad ~ log_pop + C(country)",
            data=fit,
        ).fit()
        rows.append(
            {
                "year": int(year),
                "common_loglog_slope": float(result.params["log_pop"]),
                "slope_se": float(result.bse["log_pop"]),
                "n": int(len(fit)),
                "n_countries": int(fit["CNTR_CODE"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def plot_band_comparison(bands: pd.DataFrame, interval: str, output: Path):
    fig, ax = plt.subplots(figsize=(9.6, 6.1))
    labels = list(dict.fromkeys(bands["density_band"].tolist()))
    x = np.arange(len(labels), dtype=float)

    for scenario, g in bands.groupby("scenario", sort=False):
        gg = g.set_index("density_band").reindex(labels)
        ax.plot(
            x,
            gg["pct_radiance_change_per_year"],
            marker="o",
            label=scenario.replace("_", " "),
        )

    ax.axhline(0.0, linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Reference population-density band (people km⁻²)")
    ax.set_ylabel("Estimated radiance change (% yr⁻¹)")
    ax.set_title(f"Calibration sensitivity of density-dependent growth: {interval}")
    ax.grid(alpha=0.2, axis="y")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_annual_slopes(slopes: pd.DataFrame, output: Path):
    fig, ax = plt.subplots(figsize=(9.5, 5.9))
    for scenario, g in slopes.groupby("scenario", sort=False):
        g = g.sort_values("year")
        ax.plot(
            g["year"],
            g["common_loglog_slope"],
            marker="o",
            label=scenario.replace("_", " "),
        )
    ax.axvline(2016.5, linestyle="--", linewidth=1.0)
    ax.set_xlabel("Year")
    ax.set_ylabel("Country-adjusted log10 radiance ~ log10 density slope")
    ax.set_title("Does population-radiance slope rotation survive global calibration?")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["global_calibration_rotation_test"]
    processed = configured_path(config, "processed")

    offsets, country_offsets = estimate_global_offsets(config)
    merged = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()
    data = add_scenarios(merged, offsets, config)

    scenario_cols = {
        "raw_uncorrected": "rad_raw_uncorrected",
        "legacy_five_country": "rad_legacy_five_country",
        "published_global_0p15": "rad_published_global_0p15",
        "empirical_global_dark20": "rad_empirical_global_dark20",
        "empirical_country_median_dark20":
            "rad_empirical_country_median_dark20",
    }

    # Common positive support makes scenario comparisons depend on the
    # calibration correction, not on different log-transform sample sizes.
    common_positive = np.ones(len(data), dtype=bool)
    for col in scenario_cols.values():
        common_positive &= np.isfinite(data[col]) & (data[col] > 0)
    common_positive = pd.Series(common_positive, index=data.index)

    band_frames = []
    test_rows = []
    reference_frames = []
    annual_slope_frames = []
    support_rows = []

    intervals = [
        (
            f"{cfg['first_year']}-{cfg['last_year']}",
            int(cfg["first_year"]),
            int(cfg["last_year"]),
            None,
        ),
        (
            f"{cfg['fixed_geography_first_year']}-"
            f"{cfg['fixed_geography_last_year']}_NUTS"
            f"{cfg['fixed_geography_nuts_release']}",
            int(cfg["fixed_geography_first_year"]),
            int(cfg["fixed_geography_last_year"]),
            int(cfg["fixed_geography_nuts_release"]),
        ),
    ]

    for interval_name, first_year, last_year, release in intervals:
        for scenario, rad_col in scenario_cols.items():
            raw_subset = data[
                data["year"].between(first_year, last_year, inclusive="both")
            ].copy()
            if release is not None:
                raw_subset = raw_subset[
                    raw_subset["nuts_release"] == release
                ].copy()

            support_rows.append(
                {
                    "interval": interval_name,
                    "scenario": scenario,
                    "n_rows_before_positive_filter": len(raw_subset),
                    "n_positive_scenario": int(
                        (np.isfinite(raw_subset[rad_col]) & (raw_subset[rad_col] > 0)).sum()
                    ),
                    "n_common_positive_all_scenarios": int(
                        common_positive.reindex(raw_subset.index).fillna(False).sum()
                    ),
                }
            )

            d = prepare_scenario(
                data,
                rad_col,
                config,
                first_year,
                last_year,
                release,
                common_positive_mask=common_positive,
            )

            bands, binfo = fit_band_model(
                d,
                cfg["density_band_labels"],
            )
            bands.insert(0, "interval", interval_name)
            bands.insert(1, "scenario", scenario)
            band_frames.append(bands)

            refs, cinfo = fit_continuous_model(
                d,
                cfg["reference_population_densities"],
            )
            refs.insert(0, "interval", interval_name)
            refs.insert(1, "scenario", scenario)
            reference_frames.append(refs)

            test_rows.append(
                {
                    "interval": interval_name,
                    "scenario": scenario,
                    "n_obs": len(d),
                    "n_nuts3": d["NUTS_ID"].nunique(),
                    "n_countries": d["CNTR_CODE"].nunique(),
                    **binfo,
                    **cinfo,
                }
            )

            slopes = annual_common_slopes(d)
            slopes.insert(0, "interval", interval_name)
            slopes.insert(1, "scenario", scenario)
            annual_slope_frames.append(slopes)

    bands = pd.concat(band_frames, ignore_index=True)
    tests = pd.DataFrame(test_rows)
    references = pd.concat(reference_frames, ignore_index=True)
    annual_slopes = pd.concat(annual_slope_frames, ignore_index=True)
    support = pd.DataFrame(support_rows)

    offsets_path = processed / "global_calibration_offset_scenarios.csv"
    country_offsets_path = processed / "global_calibration_country_diagnostics.csv"
    bands_path = processed / "global_calibration_density_band_growth.csv"
    tests_path = processed / "global_calibration_rotation_tests.csv"
    refs_path = processed / "global_calibration_reference_density_growth.csv"
    slopes_path = processed / "global_calibration_annual_common_slopes.csv"
    support_path = processed / "global_calibration_common_support.csv"
    summary_path = processed / "global_calibration_rotation_summary.txt"

    offsets.to_csv(offsets_path, index=False)
    country_offsets.to_csv(country_offsets_path, index=False)
    bands.to_csv(bands_path, index=False)
    tests.to_csv(tests_path, index=False)
    references.to_csv(refs_path, index=False)
    annual_slopes.to_csv(slopes_path, index=False)
    support.to_csv(support_path, index=False)

    outdir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "global_calibration_rotation"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    full_label = f"{cfg['first_year']}-{cfg['last_year']}"
    fixed_label = (
        f"{cfg['fixed_geography_first_year']}-"
        f"{cfg['fixed_geography_last_year']}_NUTS"
        f"{cfg['fixed_geography_nuts_release']}"
    )
    full_plot = outdir / "density_band_growth_global_calibration_full.pdf"
    fixed_plot = outdir / "density_band_growth_global_calibration_fixed_nuts.pdf"
    slope_plot = outdir / "annual_population_radiance_slope_calibration_comparison.pdf"

    plot_band_comparison(
        bands[bands["interval"] == full_label],
        full_label,
        full_plot,
    )
    plot_band_comparison(
        bands[bands["interval"] == fixed_label],
        fixed_label,
        fixed_plot,
    )
    plot_annual_slopes(
        annual_slopes[annual_slopes["interval"] == full_label],
        slope_plot,
    )

    with summary_path.open("w", encoding="utf-8") as h:
        h.write("GLOBAL 2017 VIIRS CALIBRATION: DOES THE DENSITY ROTATION SURVIVE?\\n")
        h.write("=" * 78 + "\\n\\n")
        h.write(
            "Diagnostic only: the production calibrated table is NOT "
            "overwritten by this stage.\\n\\n"
        )
        h.write("GLOBAL OFFSET SCENARIOS\\n")
        h.write(
            offsets.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\\n\\n")

        for interval in [full_label, fixed_label]:
            h.write(f"INTERVAL: {interval}\\n")
            h.write("-" * 72 + "\\n")
            h.write("DENSITY-BAND GROWTH (% yr-1)\\n")
            pivot = (
                bands[bands["interval"] == interval]
                .pivot(
                    index="density_band",
                    columns="scenario",
                    values="pct_radiance_change_per_year",
                )
            )
            h.write(
                pivot.to_string(
                    float_format=lambda x: f"{x:.6g}"
                )
            )
            h.write("\\n\\nROTATION TESTS\\n")
            h.write(
                tests[tests["interval"] == interval][
                    [
                        "scenario",
                        "n_obs",
                        "band_year_lr",
                        "band_year_df",
                        "band_year_p",
                        "continuous_density_year_lr",
                        "continuous_density_year_df",
                        "continuous_density_year_p",
                    ]
                ].to_string(
                    index=False,
                    float_format=lambda x: f"{x:.6g}",
                )
            )
            h.write("\\n\\nREFERENCE-DENSITY GROWTH (% yr-1)\\n")
            pivot2 = (
                references[references["interval"] == interval]
                .pivot(
                    index="population_density",
                    columns="scenario",
                    values="pct_radiance_change_per_year",
                )
            )
            h.write(
                pivot2.to_string(
                    float_format=lambda x: f"{x:.6g}"
                )
            )
            h.write("\\n\\n")

        h.write("INTERPRETATION\\n")
        h.write(
            "- The key question is whether low-density growth remains larger "
            "than high-density growth under BOTH global correction scenarios.\\n"
        )
        h.write(
            "- Strong band x year and continuous density x year tests under "
            "global correction mean the population-radiance relationship still "
            "rotates through time.\\n"
        )
        h.write(
            "- Agreement between the published 0.15 correction and the two "
            "empirical Europe-wide offsets would show that the conclusion is "
            "not sensitive to the exact global additive value.\\n"
        )
        h.write(
            "- The fixed NUTS2016 2013-2020 analysis is the cleaner geography "
            "sensitivity test.\\n"
        )
        h.write(
            "- Because Annual VNL V2 average-masked rasters contain a lit-mask "
            "operation, this NUTS-mean additive test is a robustness diagnostic. "
            "Do not promote a production correction until its behaviour has "
            "been checked against the product masking semantics.\\n"
        )

    print("\\nGLOBAL OFFSET SCENARIOS")
    print(offsets.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

    for interval in [full_label, fixed_label]:
        print(f"\\n{interval}: DENSITY-BAND GROWTH (% yr-1)")
        print(
            bands[bands["interval"] == interval]
            .pivot(
                index="density_band",
                columns="scenario",
                values="pct_radiance_change_per_year",
            )
            .to_string(float_format=lambda x: f"{x:.6g}")
        )
        print(f"\\n{interval}: ROTATION TESTS")
        print(
            tests[tests["interval"] == interval][
                [
                    "scenario",
                    "band_year_lr",
                    "band_year_df",
                    "band_year_p",
                    "continuous_density_year_lr",
                    "continuous_density_year_df",
                    "continuous_density_year_p",
                ]
            ].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )

    print("\\nWROTE")
    for p in [
        offsets_path,
        country_offsets_path,
        bands_path,
        tests_path,
        refs_path,
        slopes_path,
        support_path,
        summary_path,
        full_plot,
        fixed_plot,
        slope_plot,
    ]:
        print(p)


if __name__ == "__main__":
    main()
