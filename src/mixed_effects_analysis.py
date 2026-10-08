#!/usr/bin/env python3
"""European NUTS3 longitudinal linear mixed-effects analysis.

The response is continuous log10 corrected VIIRS radiance, so this is a
Gaussian linear mixed-effects model (LMM), i.e. the Gaussian special case of
the broader GLMM framework.

Production full model:
    log10_radiance ~ log10_population_density
                   + year_centered
                   + C(country)
                   + log10_population_density:C(country)
                   + year_centered:C(country)
                   + (1 | NUTS_ID)

NUTS3 is a random intercept. Sequential ML fits quantify added contributions
of year, country, country-specific population scaling and country-specific
temporal trends.
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


def fit_model(formula: str, data: pd.DataFrame):
    model = smf.mixedlm(
        formula,
        data=data,
        groups=data["NUTS_ID"],
        re_formula="1",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("always")
        result = model.fit(
            reml=False,
            method="lbfgs",
            maxiter=1000,
            disp=False,
        )
    return result


def nakagawa_r2(result) -> tuple[float, float, float, float, float]:
    """Nakagawa-style marginal and conditional R2 for random-intercept LMM."""
    fixed_pred = np.asarray(result.model.exog @ result.fe_params, dtype=float)
    var_fixed = float(np.var(fixed_pred, ddof=1))
    var_random = float(np.asarray(result.cov_re)[0, 0])
    var_resid = float(result.scale)
    total = var_fixed + var_random + var_resid

    if total <= 0:
        return np.nan, np.nan, var_fixed, var_random, var_resid

    marginal = var_fixed / total
    conditional = (var_fixed + var_random) / total
    return marginal, conditional, var_fixed, var_random, var_resid


def model_row(name, formula, result):
    marginal_r2, conditional_r2, var_fixed, var_random, var_resid = (
        nakagawa_r2(result)
    )
    return {
        "model": name,
        "formula": formula,
        "n_obs": int(result.nobs),
        "n_fixed_parameters": len(result.fe_params),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "bic": float(result.bic),
        "marginal_r2": marginal_r2,
        "conditional_r2": conditional_r2,
        "fixed_effect_variance": var_fixed,
        "random_intercept_variance": var_random,
        "residual_variance": var_resid,
        "converged": bool(result.converged),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument(
        "--min-country-regions",
        type=int,
        default=3,
        help="Countries with fewer unique NUTS3 regions are excluded.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")

    data = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()

    data = data[
        np.isfinite(data["log10_radiance_mean_corrected"])
        & np.isfinite(data["log10_population_density"])
    ].copy()

    counts = (
        data.groupby("CNTR_CODE")["NUTS_ID"]
        .nunique()
        .rename("n_regions")
    )
    keep_countries = counts[
        counts >= args.min_country_regions
    ].index.tolist()

    data = data[data["CNTR_CODE"].isin(keep_countries)].copy()
    data["country"] = data["CNTR_CODE"].astype("category")
    data["year_centered"] = data["year"] - data["year"].mean()
    data["log_pop"] = data["log10_population_density"]
    data["log_rad"] = data["log10_radiance_mean_corrected"]

    # Standardized continuous predictors are retained for effect-size reporting.
    data["log_pop_z"] = (
        data["log_pop"] - data["log_pop"].mean()
    ) / data["log_pop"].std(ddof=1)
    data["year_z"] = (
        data["year_centered"] - data["year_centered"].mean()
    ) / data["year_centered"].std(ddof=1)

    # Require a NUTS_ID to have at least two observations to inform the
    # longitudinal random intercept.
    n_per_region = data.groupby("NUTS_ID").size()
    repeated_ids = n_per_region[n_per_region >= 2].index
    data = data[data["NUTS_ID"].isin(repeated_ids)].copy()

    formulas = [
        (
            "M0_population",
            "log_rad ~ log_pop",
        ),
        (
            "M1_plus_year",
            "log_rad ~ log_pop + year_centered",
        ),
        (
            "M2_plus_country",
            "log_rad ~ log_pop + year_centered + C(country)",
        ),
        (
            "M3_country_density_slopes",
            "log_rad ~ log_pop + year_centered + C(country) + log_pop:C(country)",
        ),
        (
            "M4_country_year",
            "log_rad ~ log_pop + year_centered + C(country) + "
            "log_pop:C(country) + year_centered:C(country)",
        ),
        (
            "M5_plus_population_year",
            "log_rad ~ log_pop + year_centered + C(country) + "
            "log_pop:C(country) + year_centered:C(country) + "
            "log_pop:year_centered",
        ),
        (
            "M6_full_three_way",
            "log_rad ~ log_pop * year_centered * C(country)",
        ),
    ]

    results = []
    comparison_rows = []
    previous = None
    previous_name = None

    print(
        f"Rows used: {len(data):,}; countries: {data['country'].nunique()}; "
        f"repeated NUTS3 groups: {data['NUTS_ID'].nunique():,}"
    )

    for name, formula in formulas:
        print(f"\nFitting {name}: {formula}", flush=True)
        result = fit_model(formula, data)
        results.append((name, formula, result))

        if previous is not None:
            lr = 2 * (result.llf - previous.llf)
            df = len(result.fe_params) - len(previous.fe_params)
            p = chi2.sf(max(lr, 0), df) if df > 0 else np.nan
            r2_reduced = nakagawa_r2(previous)[0]
            r2_full = nakagawa_r2(result)[0]
            delta_r2 = r2_full - r2_reduced
            cohens_f2 = (
                delta_r2 / (1.0 - r2_full)
                if np.isfinite(r2_full) and r2_full < 1.0
                else np.nan
            )
            comparison_rows.append(
                {
                    "reduced_model": previous_name,
                    "full_model": name,
                    "lr_statistic": float(lr),
                    "df_added": int(df),
                    "p_value": float(p),
                    "marginal_r2_reduced": r2_reduced,
                    "marginal_r2_full": r2_full,
                    "delta_marginal_r2": delta_r2,
                    "cohens_f2_incremental": cohens_f2,
                }
            )

        previous = result
        previous_name = name

    model_table = pd.DataFrame(
        [model_row(name, formula, result) for name, formula, result in results]
    )
    comparison = pd.DataFrame(comparison_rows)

    final_name, final_formula, final = results[-1]

    standardized_formula = (
        "log_rad ~ log_pop_z * year_z * C(country)"
    )
    print(
        f"\nFitting standardized full model for comparable continuous "
        f"effect sizes: {standardized_formula}",
        flush=True,
    )
    standardized = fit_model(standardized_formula, data)

    fixed = pd.DataFrame(
        {
            "term": final.fe_params.index,
            "estimate": final.fe_params.values,
            "std_error": final.bse_fe.values,
            "z": final.fe_params.values / final.bse_fe.values,
        }
    )
    fixed["p_value"] = 2 * norm.sf(np.abs(fixed["z"]))

    random_var = float(np.asarray(final.cov_re)[0, 0])
    residual_var = float(final.scale)
    icc = random_var / (random_var + residual_var)

    standardized_fixed = pd.DataFrame(
        {
            "term": standardized.fe_params.index,
            "standardized_estimate": standardized.fe_params.values,
            "std_error": standardized.bse_fe.values,
            "z": standardized.fe_params.values / standardized.bse_fe.values,
        }
    )
    standardized_fixed["p_value"] = 2 * norm.sf(
        np.abs(standardized_fixed["z"])
    )

    variance = pd.DataFrame(
        [
            {
                "random_effect": "NUTS_ID_intercept",
                "variance": random_var,
                "residual_variance": residual_var,
                "icc_nuts3": icc,
                "n_nuts3_groups": data["NUTS_ID"].nunique(),
                "n_countries": data["country"].nunique(),
                "n_observations": len(data),
            }
        ]
    )

    model_path = processed / "mixed_model_comparison.csv"
    compare_path = processed / "mixed_model_likelihood_ratio_tests.csv"
    fixed_path = processed / "mixed_model_fixed_effects.csv"
    variance_path = processed / "mixed_model_variance_components.csv"
    standardized_path = processed / "mixed_model_standardized_fixed_effects.csv"
    summary_path = processed / "mixed_model_full_summary.txt"

    model_table.to_csv(model_path, index=False)
    comparison.to_csv(compare_path, index=False)
    fixed.to_csv(fixed_path, index=False)
    variance.to_csv(variance_path, index=False)
    standardized_fixed.to_csv(standardized_path, index=False)

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(final.summary().as_text())
        handle.write("\n\nSTANDARDIZED FULL MODEL\n")
        handle.write(standardized.summary().as_text())
        handle.write("\n\n")
        handle.write(f"NUTS3 random-intercept variance: {random_var:.8f}\n")
        handle.write(f"Residual variance: {residual_var:.8f}\n")
        handle.write(f"NUTS3 ICC: {icc:.8f}\n")

    print("\nMODEL COMPARISON")
    print(model_table.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

    print("\nSEQUENTIAL LIKELIHOOD-RATIO TESTS")
    print(comparison.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

    print("\nKEY EFFECT TESTS")
    labels = {
        "M1_plus_year": "Year",
        "M2_plus_country": "Country",
        "M3_country_density_slopes": "Population density × country",
        "M4_country_year": "Year × country",
        "M5_plus_population_year": "Population density × year",
        "M6_full_three_way": "Population density × year × country",
    }
    key = comparison.copy()
    key["effect_added"] = key["full_model"].map(labels)
    print(
        key[
            [
                "effect_added",
                "lr_statistic",
                "df_added",
                "p_value",
                "delta_marginal_r2",
                "cohens_f2_incremental",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.6g}")
    )

    print("\nNUTS3 RANDOM EFFECT")
    print(variance.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

    print()
    print(f"Wrote: {model_path}")
    print(f"Wrote: {compare_path}")
    print(f"Wrote: {fixed_path}")
    print(f"Wrote: {variance_path}")
    print(f"Wrote: {standardized_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
