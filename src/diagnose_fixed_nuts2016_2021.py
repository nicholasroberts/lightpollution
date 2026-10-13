#!/usr/bin/env python3
"""Diagnose the apparent 2021 discontinuity in the fixed-NUTS2016 panel.

This stage deliberately reuses the validated fixed-geography outputs from
validate_fixed_nuts2016_2014_2024.py. It does NOT re-extract VIIRS rasters and
does NOT define a new cohort.

Scientific purpose
------------------
The fixed NUTS2016 2014-2024 validation showed that extending the same regions
from 2020 to 2024 materially changes the fitted density-band trends. The next
question is therefore temporal: is the apparent change driven by an unusual
2020->2021 within-region pulse, by a persistent post-2021 change in slope, by
country-specific behaviour, or by density-dependent behaviour?

Working scale
-------------
Radiance is analysed as log10(calibrated VIIRS DNB radiance). First
differences remove each NUTS3 region's time-invariant level:

    delta_it = log10(R_it) - log10(R_i,t-1)

On this scale:
- a constant linear trend in levels is a constant annual delta;
- a 2021 level step is a one-year pulse for the transition ending in 2021;
- a post-2021 slope change is a persistent change in annual delta from 2022 on.

Because the response is first-differenced, the time-invariant NUTS3 intercept
has already been removed algebraically. Re-fitting a NUTS_ID random intercept
to these differences is therefore unnecessary and can put its variance on the
boundary. The diagnostic models use OLS with country fixed effects and
NUTS_ID-cluster-robust covariance, which accounts for repeated annual changes
within each region without re-introducing a redundant region intercept.

Pre-specified coefficient blocks are tested with cluster-robust Wald tests.
AIC/BIC are retained only as descriptive likelihood summaries, not as the
basis for inference.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import prepare_analysis_panel


def load_panel(config: dict) -> pd.DataFrame:
    """Load the existing fixed-NUTS2016 extraction and rebuild the same panel."""
    cfg = config["analysis"]["fixed_nuts2016_2014_2024"]
    diag = config["analysis"]["fixed_nuts2016_2021_diagnostic"]
    processed = configured_path(config, "processed")

    fixed_path = processed / diag["fixed_radiance_file"]
    reference_path = processed / diag["population_reference_file"]

    if not fixed_path.exists():
        raise FileNotFoundError(
            f"Missing {fixed_path}. Run: "
            "python src/run_pipeline.py future-radiance-fixed-2014-2024"
        )
    if not reference_path.exists():
        raise FileNotFoundError(
            f"Missing {reference_path}. Run: "
            "python src/run_pipeline.py future-radiance-fixed-2014-2024"
        )

    fixed = pd.read_parquet(fixed_path).copy()
    reference = pd.read_csv(reference_path).copy()

    panel = prepare_analysis_panel(fixed, reference, config)

    # The population reference is stored as CSV, so pandas reads density_band
    # back as plain strings and loses the ordered categorical metadata created
    # by pd.cut(). Restore the configured scientific order explicitly. Without
    # this, Patsy chooses the alphabetically first label ("100-300") as the
    # model reference category, while downstream contrasts assume "<30".
    density_labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )
    panel["density_band"] = pd.Categorical(
        panel["density_band"],
        categories=density_labels,
        ordered=True,
    )
    if panel["density_band"].isna().any():
        raise RuntimeError(
            "At least one density-band value is not in the configured "
            f"ordered labels: {density_labels}"
        )

    panel = panel.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)

    first_year = int(cfg["first_year"])
    last_year = int(cfg["last_year"])
    expected_years = list(range(first_year, last_year + 1))

    observed_years = sorted(panel["year"].astype(int).unique().tolist())
    if observed_years != expected_years:
        raise RuntimeError(
            f"Fixed panel years are {observed_years}; expected {expected_years}."
        )

    n_per_region = panel.groupby("NUTS_ID")["year"].nunique()
    if not (n_per_region == len(expected_years)).all():
        raise RuntimeError(
            "2021 diagnostic requires the same complete fixed-NUTS2016 panel "
            "used by the validation stage."
        )

    expected_n = diag.get("expected_n_nuts3")
    if expected_n is not None and panel["NUTS_ID"].nunique() != int(expected_n):
        raise RuntimeError(
            f"Panel contains {panel['NUTS_ID'].nunique()} NUTS3 regions; "
            f"configured expectation is {int(expected_n)}. "
            "Do not continue until the cohort mismatch is explained."
        )

    return panel


def make_first_differences(panel: pd.DataFrame) -> pd.DataFrame:
    """Create within-region annual changes on the log10 radiance scale."""
    d = panel[
        [
            "NUTS_ID",
            "CNTR_CODE",
            "density_band",
            "reference_population_density",
            "year",
            "radiance_mean_corrected",
            "log_rad",
        ]
    ].copy()
    d = d.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)

    d["previous_year"] = d.groupby("NUTS_ID")["year"].shift(1)
    d["previous_log_rad"] = d.groupby("NUTS_ID")["log_rad"].shift(1)
    d["previous_radiance"] = d.groupby("NUTS_ID")[
        "radiance_mean_corrected"
    ].shift(1)

    d["delta_log10"] = d["log_rad"] - d["previous_log_rad"]
    d["pct_change"] = 100.0 * (10.0 ** d["delta_log10"] - 1.0)
    d = d[d["previous_year"].notna()].copy()

    bad_gap = d["year"].astype(int) - d["previous_year"].astype(int) != 1
    if bad_gap.any():
        raise RuntimeError("Non-consecutive years found in complete panel.")

    d["end_year"] = d["year"].astype(int)
    d["jump_2021"] = (d["end_year"] == 2021).astype(int)
    d["post_2021"] = (d["end_year"] >= 2022).astype(int)
    d["country"] = d["CNTR_CODE"].astype("category")
    # Preserve the already-restored ordered density categories from load_panel.
    if not isinstance(d["density_band"].dtype, pd.CategoricalDtype):
        raise RuntimeError(
            "density_band lost its configured categorical ordering before "
            "model fitting."
        )
    return d


def summarise_transition(g: pd.DataFrame) -> dict:
    x = g["delta_log10"].to_numpy(float)
    pct = g["pct_change"].to_numpy(float)
    return {
        "n_nuts3": int(g["NUTS_ID"].nunique()),
        "n_countries": int(g["CNTR_CODE"].nunique()),
        "mean_delta_log10": float(np.mean(x)),
        "median_delta_log10": float(np.median(x)),
        "se_delta_log10": (
            float(np.std(x, ddof=1) / np.sqrt(len(x)))
            if len(x) > 1
            else np.nan
        ),
        "geometric_pct_change": float(100.0 * (10.0 ** np.mean(x) - 1.0)),
        "median_region_pct_change": float(np.median(pct)),
        "pct_regions_increasing": float(100.0 * np.mean(x > 0)),
    }


def annual_summary(diff: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for year, g in diff.groupby("end_year", observed=True):
        row = {"end_year": int(year), **summarise_transition(g)}
        rows.append(row)
    return pd.DataFrame(rows).sort_values("end_year")


def density_summary(diff: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, band), g in diff.groupby(
        ["end_year", "density_band"], observed=True
    ):
        row = {
            "end_year": int(year),
            "density_band": str(band),
            **summarise_transition(g),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["end_year", "density_band"])


def country_summary(diff: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (year, country), g in diff.groupby(
        ["end_year", "CNTR_CODE"], observed=True
    ):
        row = {
            "end_year": int(year),
            "CNTR_CODE": str(country),
            **summarise_transition(g),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["end_year", "CNTR_CODE"])


def country_balanced_summary(country: pd.DataFrame) -> pd.DataFrame:
    """Summarise country medians equally, so large countries do not dominate."""
    rows = []
    for year, g in country.groupby("end_year"):
        vals = g["median_delta_log10"].to_numpy(float)
        rows.append(
            {
                "end_year": int(year),
                "n_countries": int(len(g)),
                "mean_country_median_delta_log10": float(np.mean(vals)),
                "median_country_median_delta_log10": float(np.median(vals)),
                "country_balanced_geometric_pct_change": float(
                    100.0 * (10.0 ** np.mean(vals) - 1.0)
                ),
                "pct_countries_positive": float(100.0 * np.mean(vals > 0)),
            }
        )
    return pd.DataFrame(rows).sort_values("end_year")


def fit_cluster_ols(formula: str, data: pd.DataFrame):
    """OLS point estimates with covariance clustered by fixed NUTS3 region."""
    model = smf.ols(formula, data=data)
    return model.fit(
        cov_type="cluster",
        cov_kwds={"groups": data["NUTS_ID"], "use_correction": True},
    )


def robust_wald(result, restrictions: list[dict[str, float]]):
    """Joint cluster-robust Wald chi-square test of linear restrictions."""
    names = list(result.params.index)
    rows = []
    for weights in restrictions:
        row = np.zeros(len(names), dtype=float)
        for term, value in weights.items():
            if term not in names:
                raise RuntimeError(
                    f"Restriction term {term!r} not found in fitted model."
                )
            row[names.index(term)] = float(value)
        rows.append(row)
    R = np.vstack(rows)
    test = result.wald_test(R, use_f=False, scalar=True)
    return float(test.statistic), int(R.shape[0]), float(test.pvalue)


def cluster_contrast(result, weights: dict[str, float]):
    """Estimate one linear contrast using the cluster-robust covariance."""
    names = list(result.params.index)
    L = np.zeros(len(names), dtype=float)
    for term, value in weights.items():
        if term not in names:
            raise RuntimeError(f"Contrast term {term!r} not found in model.")
        L[names.index(term)] = float(value)
    beta = result.params.to_numpy(float)
    cov = result.cov_params().loc[names, names].to_numpy(float)
    est = float(L @ beta)
    se = float(np.sqrt(max(0.0, L @ cov @ L)))
    z = est / se if se > 0 else np.nan
    from scipy.stats import norm
    p = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return est, se, p


def _interaction_terms(result, variable: str, factor: str) -> list[str]:
    names = list(result.params.index)
    return [
        name
        for name in names
        if variable in name and f"C({factor})" in name and ":" in name
    ]


def _categorical_piecewise_restrictions(result, diff: pd.DataFrame):
    """Restrictions implied by a constant/pulse/post-2021 delta model."""
    years = sorted(diff["end_year"].astype(int).unique().tolist())
    reference = years[0]
    if reference != 2015 or 2021 not in years or 2022 not in years:
        raise RuntimeError(
            f"Unexpected transition years for 2021 diagnostic: {years}"
        )

    restrictions = []
    # Before the 2021 pulse, every transition must equal the 2015 reference.
    for year in [y for y in years if reference < y < 2021]:
        restrictions.append({f"C(end_year)[T.{year}]": 1.0})

    # Post-2021 transitions are all constrained to the same shifted mean.
    post_reference = "C(end_year)[T.2022]"
    for year in [y for y in years if y > 2022]:
        restrictions.append(
            {
                f"C(end_year)[T.{year}]": 1.0,
                post_reference: -1.0,
            }
        )
    return restrictions


def fit_diagnostic_models(diff: pd.DataFrame):
    """Fit the pre-specified first-difference models with clustered covariance."""
    base_formula = "delta_log10 ~ C(country) + C(density_band)"
    pulse_formula = base_formula + " + jump_2021"
    piecewise_formula = pulse_formula + " + post_2021"
    density_formula = (
        piecewise_formula
        + " + jump_2021:C(density_band)"
        + " + post_2021:C(density_band)"
    )
    country_formula = (
        density_formula
        + " + jump_2021:C(country)"
        + " + post_2021:C(country)"
    )
    categorical_formula = (
        "delta_log10 ~ C(country) + C(density_band) + C(end_year)"
    )

    formulas = {
        "constant_change": base_formula,
        "plus_2021_pulse": pulse_formula,
        "plus_post_2021_change": piecewise_formula,
        "plus_density_specific_breaks": density_formula,
        "plus_country_specific_breaks": country_formula,
        "categorical_transition_year": categorical_formula,
    }
    models = {
        name: fit_cluster_ols(formula, diff)
        for name, formula in formulas.items()
    }

    test_specs = [
        (
            "2021 pulse beyond constant annual change",
            "constant_change",
            "plus_2021_pulse",
            [{"jump_2021": 1.0}],
        ),
        (
            "persistent post-2021 change beyond one-year pulse",
            "plus_2021_pulse",
            "plus_post_2021_change",
            [{"post_2021": 1.0}],
        ),
    ]

    density_model = models["plus_density_specific_breaks"]
    density_terms = (
        _interaction_terms(density_model, "jump_2021", "density_band")
        + _interaction_terms(density_model, "post_2021", "density_band")
    )
    test_specs.append(
        (
            "density-specific 2021/post-2021 behaviour",
            "plus_post_2021_change",
            "plus_density_specific_breaks",
            [{term: 1.0} for term in density_terms],
        )
    )

    country_model = models["plus_country_specific_breaks"]
    country_terms = (
        _interaction_terms(country_model, "jump_2021", "country")
        + _interaction_terms(country_model, "post_2021", "country")
    )
    test_specs.append(
        (
            "country-specific 2021/post-2021 behaviour",
            "plus_density_specific_breaks",
            "plus_country_specific_breaks",
            [{term: 1.0} for term in country_terms],
        )
    )

    categorical_model = models["categorical_transition_year"]
    test_specs.append(
        (
            "remaining year-to-year structure beyond simple 2021 piecewise model",
            "plus_post_2021_change",
            "categorical_transition_year",
            _categorical_piecewise_restrictions(categorical_model, diff),
        )
    )

    tests = []
    for label, reduced, full, restrictions in test_specs:
        stat, df_wald, p = robust_wald(models[full], restrictions)
        tests.append(
            {
                "test": label,
                "reduced_model": reduced,
                "full_model": full,
                "wald_chi2": stat,
                "df": df_wald,
                "p": p,
                "reduced_aic": float(models[reduced].aic),
                "full_aic": float(models[full].aic),
            }
        )

    model_info = []
    n_clusters = int(diff["NUTS_ID"].nunique())
    for name, model in models.items():
        model_info.append(
            {
                "model": name,
                "estimator": "OLS with NUTS_ID-cluster-robust covariance",
                "n_obs": int(model.nobs),
                "n_clusters": n_clusters,
                "n_fixed_effects": int(len(model.params)),
                "log_likelihood": float(model.llf),
                "aic": float(model.aic),
                "bic": float(model.bic),
                "covariance": str(model.cov_type),
            }
        )

    return models, pd.DataFrame(tests), pd.DataFrame(model_info)


def coefficient_table(models: dict) -> pd.DataFrame:
    rows = []
    keep_models = [
        "plus_post_2021_change",
        "plus_density_specific_breaks",
        "plus_country_specific_breaks",
        "categorical_transition_year",
    ]
    for model_name in keep_models:
        model = models[model_name]
        fe = model.params
        se = model.bse
        for term in fe.index:
            if (
                term == "Intercept"
                or "jump_2021" in term
                or "post_2021" in term
                or "C(end_year)" in term
            ):
                est = float(fe[term])
                s = float(se[term])
                rows.append(
                    {
                        "model": model_name,
                        "term": term,
                        "estimate_log10": est,
                        "se": s,
                        "ci_low": est - 1.96 * s,
                        "ci_high": est + 1.96 * s,
                    }
                )
    return pd.DataFrame(rows)


def density_break_contrasts(model, labels: list[str]) -> pd.DataFrame:
    """Within-country density contrasts for the 2021 pulse and post shift."""
    rows = []
    ref = labels[0]
    for band in labels:
        jump_weights = {"jump_2021": 1.0}
        post_weights = {"post_2021": 1.0}
        if band != ref:
            for base, weights in [
                ("jump_2021", jump_weights),
                ("post_2021", post_weights),
            ]:
                candidates = [
                    f"{base}:C(density_band)[T.{band}]",
                    f"C(density_band)[T.{band}]:{base}",
                ]
                found = next(
                    (x for x in candidates if x in model.params.index),
                    None,
                )
                if found is None:
                    raise RuntimeError(
                        f"Could not find density interaction for {band}: {base}"
                    )
                weights[found] = 1.0

        jump, jump_se, jump_p = cluster_contrast(model, jump_weights)
        post, post_se, post_p = cluster_contrast(model, post_weights)
        rows.append(
            {
                "density_band": band,
                "reference_band": ref,
                "model": "plus_density_specific_breaks",
                "estimated_2021_pulse_log10": jump,
                "pulse_se": jump_se,
                "pulse_p": jump_p,
                "estimated_post_2021_delta_shift_log10": post,
                "post_shift_se": post_se,
                "post_shift_p": post_p,
            }
        )
    return pd.DataFrame(rows)


def make_plots(
    annual: pd.DataFrame,
    density: pd.DataFrame,
    country: pd.DataFrame,
    labels: list[str],
    output_dir: Path,
    dpi: int,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.plot(
        annual["end_year"],
        annual["geometric_pct_change"],
        marker="o",
        linewidth=1.6,
    )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(2021, linestyle="--", linewidth=1.0)
    ax.set_xlabel("Transition ending in year")
    ax.set_ylabel("Mean within-region radiance change (%)")
    ax.set_title("Annual within-region change on fixed NUTS2016 geography")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        output_dir / "annual_within_region_change.pdf",
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    for band in labels:
        g = density[density["density_band"] == band].sort_values("end_year")
        ax.plot(
            g["end_year"],
            g["geometric_pct_change"],
            marker="o",
            linewidth=1.3,
            label=band,
        )
    ax.axhline(0, linewidth=0.8)
    ax.axvline(2021, linestyle="--", linewidth=1.0)
    ax.set_xlabel("Transition ending in year")
    ax.set_ylabel("Mean within-region radiance change (%)")
    ax.set_title("Annual fixed-NUTS2016 change by fixed density band")
    ax.legend(title="People km$^{-2}$")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        output_dir / "annual_change_by_density_band.pdf",
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(fig)

    c21 = country[country["end_year"] == 2021].copy()
    c21 = c21.sort_values("geometric_pct_change")
    fig_height = max(6.0, 0.25 * len(c21) + 1.5)
    fig, ax = plt.subplots(figsize=(8.5, fig_height))
    y = np.arange(len(c21))
    ax.barh(y, c21["geometric_pct_change"])
    ax.axvline(0, linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(c21["CNTR_CODE"])
    ax.set_xlabel("Mean within-region change, 2020 to 2021 (%)")
    ax.set_ylabel("Country")
    ax.set_title("2020->2021 change by country on fixed NUTS2016 geography")
    ax.grid(alpha=0.2, axis="x")
    fig.tight_layout()
    fig.savefig(
        output_dir / "country_change_2020_2021.pdf",
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    config = load_config(args.config)
    diag_cfg = config["analysis"]["fixed_nuts2016_2021_diagnostic"]
    labels = list(
        config["analysis"]["density_temporal_change"]["density_band_labels"]
    )

    panel = load_panel(config)
    diff = make_first_differences(panel)

    annual = annual_summary(diff)
    density = density_summary(diff)
    country = country_summary(diff)
    country_balanced = country_balanced_summary(country)

    models, tests, model_info = fit_diagnostic_models(diff)
    coefficients = coefficient_table(models)
    density_contrasts = density_break_contrasts(
        models["plus_density_specific_breaks"], labels
    )

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_2021_diagnostic"
    )
    processed.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    outputs = {
        "first_differences": processed
        / "fixed_nuts2016_annual_first_differences_2015_2024.csv",
        "annual": processed
        / "fixed_nuts2016_annual_change_summary_2015_2024.csv",
        "density": processed
        / "fixed_nuts2016_annual_change_by_density_2015_2024.csv",
        "country": processed
        / "fixed_nuts2016_annual_change_by_country_2015_2024.csv",
        "country_balanced": processed
        / "fixed_nuts2016_country_balanced_annual_change_2015_2024.csv",
        "tests": processed / "fixed_nuts2016_2021_model_tests.csv",
        "model_info": processed / "fixed_nuts2016_2021_model_info.csv",
        "coefficients": processed
        / "fixed_nuts2016_2021_model_coefficients.csv",
        "density_contrasts": processed
        / "fixed_nuts2016_2021_density_break_contrasts.csv",
        "summary": processed / "fixed_nuts2016_2021_diagnostic_summary.txt",
    }

    diff.to_csv(outputs["first_differences"], index=False)
    annual.to_csv(outputs["annual"], index=False)
    density.to_csv(outputs["density"], index=False)
    country.to_csv(outputs["country"], index=False)
    country_balanced.to_csv(outputs["country_balanced"], index=False)
    tests.to_csv(outputs["tests"], index=False)
    model_info.to_csv(outputs["model_info"], index=False)
    coefficients.to_csv(outputs["coefficients"], index=False)
    density_contrasts.to_csv(outputs["density_contrasts"], index=False)

    make_plots(
        annual,
        density,
        country,
        labels,
        figures,
        int(diag_cfg["figure_dpi"]),
    )

    c21 = country[country["end_year"] == 2021].copy()
    b21 = country_balanced[
        country_balanced["end_year"] == 2021
    ].iloc[0]

    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 DIAGNOSTIC OF THE 2021 DISCONTINUITY\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Panel: {panel['NUTS_ID'].nunique()} fixed NUTS3 regions; "
            f"{panel['CNTR_CODE'].nunique()} countries; "
            f"{panel['year'].min()}-{panel['year'].max()}.\n"
        )
        h.write(
            "Cohort and density classes are inherited unchanged from the "
            "fixed-NUTS2016 2014-2024 validation.\n"
        )
        h.write(
            "Response: within-region annual difference in log10 calibrated "
            "VIIRS DNB radiance.\n\n"
        )

        h.write("ANNUAL WITHIN-REGION CHANGE\n")
        h.write(
            annual.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nCOUNTRY-BALANCED ANNUAL CHANGE\n")
        h.write(
            country_balanced.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\n2020->2021 COUNTRY SUMMARY\n")
        h.write(
            c21.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nCLUSTER-ROBUST WALD TESTS\n")
        h.write(
            tests.to_string(
                index=False, float_format=lambda x: f"{x:.6g}"
            )
        )
        h.write("\n\nKEY 2021 DESCRIPTIVE DIAGNOSTIC\n")
        h.write(
            f"Country-balanced mean of country median changes: "
            f"{b21['country_balanced_geometric_pct_change']:+.4f}%.\n"
        )
        h.write(
            f"Countries with positive median 2020->2021 change: "
            f"{b21['pct_countries_positive']:.1f}%.\n"
        )
        h.write(
            "\nINTERPRETATION RULES (pre-specified)\n"
            "- A strong common 2021 pulse with many countries moving in the "
            "same direction points to a Europe-wide temporal discontinuity, "
            "not changing region membership.\n"
            "- A significant persistent post-2021 term indicates a slope "
            "change rather than only a one-year level pulse.\n"
            "- Density-specific break terms test whether the discontinuity "
            "depends systematically on the fixed population-density gradient.\n"
            "- Country-specific break terms test whether national temporal "
            "behaviour explains substantial additional variation.\n"
            "- If the categorical-year model substantially outperforms the "
            "simple 2021 piecewise model, the final analysis should not force "
            "the record into a single break-point formulation.\n"
        )

    print("\nLOCKED FIXED-NUTS2016 PANEL")
    print(
        f"{panel['NUTS_ID'].nunique()} NUTS3; "
        f"{panel['CNTR_CODE'].nunique()} countries; "
        f"{len(panel)} annual rows"
    )
    print("\nANNUAL WITHIN-REGION CHANGE")
    print(
        annual.to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\n2020 -> 2021 CHANGE BY DENSITY BAND")
    print(
        density[density["end_year"] == 2021].to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\n2020 -> 2021 COUNTRY MEDIANS")
    print(
        c21[
            [
                "CNTR_CODE",
                "n_nuts3",
                "geometric_pct_change",
                "median_region_pct_change",
                "pct_regions_increasing",
            ]
        ].to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\nCLUSTER-ROBUST WALD TESTS")
    print(
        tests[
            ["test", "wald_chi2", "df", "p", "reduced_aic", "full_aic"]
        ].to_string(
            index=False, float_format=lambda x: f"{x:.6g}"
        )
    )
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
