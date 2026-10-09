#!/usr/bin/env python3
"""Quantify how much common rural-urban temporal structure explains country effects.

This analysis compares the strength of three country-specific fixed-effect
blocks before and after replacing a simple linear population x year surface
with the common nonlinear density-dependent temporal surface identified by
analyze_density_temporal_change.py:

  1. country-specific temporal trends: year x country
  2. country-specific population slopes: log(population density) x country
  3. country-specific rotation: year x log(population density) x country

For each block, the number of country-specific parameters is held constant.
The key descriptive effect size is the proportional reduction in the block's
likelihood-ratio statistic and LR pseudo-R2 after the common nonlinear surface
is included. This is a decomposition diagnostic, not a causal proportion.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_density_temporal_change import fit_mixedlm, lrt, prepare_data
from project_config import configured_path, load_config


def nakagawa_r2(result):
    fixed_pred = np.asarray(result.model.exog @ result.fe_params, dtype=float)
    var_fixed = float(np.var(fixed_pred, ddof=1))
    var_random = float(np.asarray(result.cov_re)[0, 0])
    var_resid = float(result.scale)
    total = var_fixed + var_random + var_resid
    if total <= 0:
        return np.nan, np.nan
    return (
        var_fixed / total,
        (var_fixed + var_random) / total,
    )


def lr_pseudo_r2(lr, n):
    return float(1.0 - np.exp(-max(float(lr), 0.0) / float(n)))


def block_test(name, reduced_formula, full_formula, data):
    print(f"\n{name}")
    print(f"  reduced: {reduced_formula}")
    reduced, opt0 = fit_mixedlm(reduced_formula, data)
    print(f"  full:    {full_formula}")
    full, opt1 = fit_mixedlm(full_formula, data)
    lr, df, p = lrt(reduced, full)
    pseudo = lr_pseudo_r2(lr, len(data))
    return {
        "block": name,
        "reduced_formula": reduced_formula,
        "full_formula": full_formula,
        "lr": lr,
        "df": df,
        "p": p,
        "lr_pseudo_r2": pseudo,
        "reduced_optimizer": opt0,
        "full_optimizer": opt1,
        "reduced_marginal_r2": nakagawa_r2(reduced)[0],
        "full_marginal_r2": nakagawa_r2(full)[0],
        "reduced_aic": float(reduced.aic),
        "full_aic": float(full.aic),
        "reduced_bic": float(reduced.bic),
        "full_bic": float(full.bic),
    }


def model_summary(name, formula, data):
    print(f"\nWhole model: {name}")
    result, opt = fit_mixedlm(formula, data)
    marginal, conditional = nakagawa_r2(result)
    return {
        "model": name,
        "formula": formula,
        "n_fixed_parameters": len(result.fe_params),
        "log_likelihood": float(result.llf),
        "aic": float(result.aic),
        "bic": float(result.bic),
        "marginal_r2": marginal,
        "conditional_r2": conditional,
        "random_intercept_variance": float(np.asarray(result.cov_re)[0, 0]),
        "residual_variance": float(result.scale),
        "optimizer": opt,
        "converged": bool(result.converged),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    data = prepare_data(config)
    n = len(data)

    # Simple linear common surface used as the "before" reference.
    linear_common = (
        "log_rad ~ C(country) + year_c * log_pop_c"
    )

    # Common nonlinear rural-urban temporal surface from the previous analysis.
    surface_common = (
        "log_rad ~ C(country) + "
        "year_c * (log_pop_c + I(log_pop_c ** 2) + I(log_pop_c ** 3))"
    )

    rows = []

    # Country-specific temporal trend block.
    rows.append({
        "framework": "linear_reference",
        **block_test(
            "country_year",
            linear_common,
            linear_common + " + year_c:C(country)",
            data,
        ),
    })
    rows.append({
        "framework": "nonlinear_density_time_surface",
        **block_test(
            "country_year",
            surface_common,
            surface_common + " + year_c:C(country)",
            data,
        ),
    })

    # Country-specific density slope block.
    rows.append({
        "framework": "linear_reference",
        **block_test(
            "country_density",
            linear_common,
            linear_common + " + log_pop_c:C(country)",
            data,
        ),
    })
    rows.append({
        "framework": "nonlinear_density_time_surface",
        **block_test(
            "country_density",
            surface_common,
            surface_common + " + log_pop_c:C(country)",
            data,
        ),
    })

    # Country-specific rotation through time. Respect hierarchy in both
    # frameworks by including country x year and country x density first.
    linear_rotation_base = (
        linear_common
        + " + year_c:C(country)"
        + " + log_pop_c:C(country)"
    )
    surface_rotation_base = (
        surface_common
        + " + year_c:C(country)"
        + " + log_pop_c:C(country)"
    )
    rows.append({
        "framework": "linear_reference",
        **block_test(
            "country_rotation",
            linear_rotation_base,
            linear_rotation_base + " + year_c:log_pop_c:C(country)",
            data,
        ),
    })
    rows.append({
        "framework": "nonlinear_density_time_surface",
        **block_test(
            "country_rotation",
            surface_rotation_base,
            surface_rotation_base + " + year_c:log_pop_c:C(country)",
            data,
        ),
    })

    tests = pd.DataFrame(rows)

    # Pair before/after results and quantify the reduction in residual
    # country-specific structure.
    paired = []
    for block in ["country_year", "country_density", "country_rotation"]:
        a = tests[
            (tests["block"] == block)
            & (tests["framework"] == "linear_reference")
        ].iloc[0]
        b = tests[
            (tests["block"] == block)
            & (tests["framework"] == "nonlinear_density_time_surface")
        ].iloc[0]

        lr_reduction = (
            1.0 - b["lr"] / a["lr"]
            if a["lr"] > 0 else np.nan
        )
        pr2_reduction = (
            1.0 - b["lr_pseudo_r2"] / a["lr_pseudo_r2"]
            if a["lr_pseudo_r2"] > 0 else np.nan
        )

        paired.append({
            "block": block,
            "linear_reference_lr": a["lr"],
            "surface_residual_lr": b["lr"],
            "fraction_lr_removed_by_surface": lr_reduction,
            "percent_lr_removed_by_surface": 100.0 * lr_reduction,
            "linear_reference_lr_pseudo_r2": a["lr_pseudo_r2"],
            "surface_residual_lr_pseudo_r2": b["lr_pseudo_r2"],
            "fraction_pseudo_r2_removed_by_surface": pr2_reduction,
            "percent_pseudo_r2_removed_by_surface": 100.0 * pr2_reduction,
            "linear_reference_p": a["p"],
            "surface_residual_p": b["p"],
            "df_block": int(a["df"]),
        })

    decomposition = pd.DataFrame(paired)

    # Whole-model comparison: same country-specific terms, but either a simple
    # linear or common nonlinear density-time surface.
    country_heavy_linear = (
        linear_rotation_base + " + year_c:log_pop_c:C(country)"
    )
    country_heavy_surface = (
        surface_rotation_base + " + year_c:log_pop_c:C(country)"
    )
    models = pd.DataFrame([
        model_summary("country_heavy_linear", country_heavy_linear, data),
        model_summary("country_heavy_with_common_surface", country_heavy_surface, data),
        model_summary("common_surface_without_country_interactions", surface_common, data),
    ])

    processed = configured_path(config, "processed")
    figdir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "density_temporal_change"
    )
    figdir.mkdir(parents=True, exist_ok=True)

    tests_path = processed / "density_country_block_tests.csv"
    decomp_path = processed / "density_country_structure_decomposition.csv"
    models_path = processed / "density_country_model_comparison.csv"
    summary_path = processed / "density_country_structure_decomposition_summary.txt"

    tests.to_csv(tests_path, index=False)
    decomposition.to_csv(decomp_path, index=False)
    models.to_csv(models_path, index=False)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    x = np.arange(len(decomposition))
    vals = decomposition["percent_lr_removed_by_surface"].to_numpy()
    ax.bar(x, vals)
    ax.axhline(0, linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([
        "Country × year",
        "Country × density",
        "Country × density × year",
    ])
    ax.set_ylabel("Country-block LR removed by common surface (%)")
    ax.set_title(
        "How much country-specific structure is absorbed by\n"
        "the common nonlinear density × time surface"
    )
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    plot_path = figdir / "country_structure_explained_by_density_time_surface.pdf"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    with summary_path.open("w", encoding="utf-8") as h:
        h.write("COUNTRY STRUCTURE EXPLAINED BY COMMON DENSITY-TIME SURFACE\n")
        h.write("=" * 76 + "\n\n")
        h.write(
            f"Observations: {n}; repeated NUTS3: {data['NUTS_ID'].nunique()}; "
            f"countries: {data['CNTR_CODE'].nunique()}\n\n"
        )
        h.write(
            "The diagnostic compares the same country-specific block before "
            "and after fitting a common cubic log-population-density surface "
            "whose temporal trend also varies with density. Percent removed is "
            "1 - residual_block_LR / original_block_LR. It is an effect-size "
            "decomposition, not a causal proportion.\n\n"
        )
        h.write("BLOCK DECOMPOSITION\n")
        h.write(
            decomposition.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nWHOLE-MODEL COMPARISON\n")
        h.write(
            models[[
                "model", "n_fixed_parameters", "log_likelihood", "aic", "bic",
                "marginal_r2", "conditional_r2",
                "random_intercept_variance", "residual_variance",
            ]].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        h.write("\n\nINTERPRETATION\n")
        h.write(
            "- Positive percent removed means the common density-time surface "
            "has absorbed part of that country-specific structure.\n"
        )
        h.write(
            "- A value near 100% would mean little country-specific block "
            "remains; 0% means essentially no reduction; a negative value "
            "means the residual country block is stronger after the common "
            "surface is fitted.\n"
        )
        h.write(
            "- country_rotation is the most direct test of whether the common "
            "rural-urban process explains country differences in how the "
            "population-radiance slope changes through time.\n"
        )

    print("\nCOUNTRY-STRUCTURE DECOMPOSITION")
    print(
        decomposition.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nWHOLE-MODEL COMPARISON")
    print(
        models[[
            "model", "n_fixed_parameters", "log_likelihood", "aic", "bic",
            "marginal_r2", "conditional_r2",
        ]].to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nWROTE")
    for p in [tests_path, decomp_path, models_path, summary_path, plot_path]:
        print(p)


if __name__ == "__main__":
    main()
