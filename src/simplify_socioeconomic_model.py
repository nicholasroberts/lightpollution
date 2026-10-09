#!/usr/bin/env python3
"""Crawley-style top-down simplification of the socioeconomic NUTS3 LMM.

The workflow follows the maximum-model -> model-simplification -> minimum
adequate model logic described by Crawley. Fixed effects are simplified from
the highest-order interactions downwards, preserving marginality.

Maximal fixed-effects structure:
    log_pop * log_gdp_pc * year_centered * C(country)

Random-effects structure is held constant throughout:
    (1 | NUTS_ID)

All fixed-effect comparisons use maximum likelihood (ML), not REML, and use
likelihood-ratio tests between nested models. At each interaction order the
currently removable terms are tested; the least significant term is deleted
when p >= alpha, the model is refitted, and the process repeats. Lower-order
terms required by retained higher-order interactions are not eligible for
deletion.
"""

from __future__ import annotations

import argparse
import itertools
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


TOKENS = ("log_pop", "log_gdp_pc", "year_centered", "C(country)")
TOKEN_ORDER = {token: i for i, token in enumerate(TOKENS)}


def term_key(term: frozenset[str]) -> tuple:
    return (
        len(term),
        tuple(TOKEN_ORDER[x] for x in sorted(term, key=TOKEN_ORDER.get)),
    )


def term_string(term: frozenset[str]) -> str:
    ordered = sorted(term, key=TOKEN_ORDER.get)
    return ":".join(ordered)


def all_hierarchical_terms() -> set[frozenset[str]]:
    terms = set()
    for size in range(1, len(TOKENS) + 1):
        for combo in itertools.combinations(TOKENS, size):
            terms.add(frozenset(combo))
    return terms


def formula_from_terms(terms: set[frozenset[str]]) -> str:
    ordered = sorted(terms, key=term_key)
    rhs = " + ".join(term_string(term) for term in ordered)
    return f"log_rad ~ {rhs}"


def removable_terms(
    terms: set[frozenset[str]],
    order: int,
) -> list[frozenset[str]]:
    """Terms of this order with no retained strict superset."""
    candidates = []
    for term in terms:
        if len(term) != order:
            continue
        has_parent = any(
            term < other
            for other in terms
        )
        if not has_parent:
            candidates.append(term)
    return sorted(candidates, key=term_key)


def fit_model(formula: str, data: pd.DataFrame):
    model = smf.mixedlm(
        formula,
        data=data,
        groups=data["NUTS_ID"],
        re_formula="1",
    )

    failures = []
    for method in ("lbfgs", "bfgs", "cg"):
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = model.fit(
                    reml=False,
                    method=method,
                    maxiter=2000,
                    disp=False,
                )
            if result.converged:
                return result, method, [str(w.message) for w in caught]
            failures.append(f"{method}: did not converge")
        except Exception as exc:
            failures.append(f"{method}: {exc}")

    raise RuntimeError(
        "Mixed model failed with all optimizers: "
        + " | ".join(failures)
    )


def lrt(full, reduced) -> tuple[float, int, float]:
    lr = 2.0 * (full.llf - reduced.llf)
    df = len(full.fe_params) - len(reduced.fe_params)
    if df <= 0:
        return float(lr), int(df), np.nan
    p = chi2.sf(max(lr, 0.0), df)
    return float(lr), int(df), float(p)


def nakagawa_r2(result):
    fixed_pred = np.asarray(
        result.model.exog @ result.fe_params,
        dtype=float,
    )
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


def prepare_data(config: dict, min_country_regions: int) -> pd.DataFrame:
    processed = configured_path(config, "processed")
    path = processed / config["outputs"]["merged_socioeconomic"]
    data = pd.read_parquet(path).copy()

    metric = config["socioeconomic"]["model_gdp_metric"]
    log_metric = f"log10_{metric}"

    required = [
        "log10_radiance_mean_corrected",
        "log10_population_density",
        log_metric,
        "year",
        "CNTR_CODE",
        "NUTS_ID",
    ]
    missing = [c for c in required if c not in data.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    complete = np.ones(len(data), dtype=bool)
    for col in [
        "log10_radiance_mean_corrected",
        "log10_population_density",
        log_metric,
        "year",
    ]:
        complete &= np.isfinite(pd.to_numeric(data[col], errors="coerce"))

    data = data.loc[complete].copy()

    counts = data.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    keep = counts[counts >= min_country_regions].index
    data = data[data["CNTR_CODE"].isin(keep)].copy()

    data["country"] = data["CNTR_CODE"].astype("category")
    data["log_rad"] = data["log10_radiance_mean_corrected"]
    data["log_pop"] = data["log10_population_density"]
    data["log_gdp_pc"] = data[log_metric]
    data["year_centered"] = data["year"] - data["year"].mean()

    # Use a single complete-case dataset for every nested comparison.
    repeated = data.groupby("NUTS_ID").size()
    repeated_ids = repeated[repeated >= 2].index
    data = data[data["NUTS_ID"].isin(repeated_ids)].copy()

    # Standardize continuous predictors for conditioning only. The model is
    # expressed in log units / years for interpretability, but centering GDP
    # and population reduces correlation among interaction columns.
    for source, target in [
        ("log_pop", "log_pop"),
        ("log_gdp_pc", "log_gdp_pc"),
    ]:
        mean = data[source].mean()
        data[target] = data[source] - mean

    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--alpha", type=float, default=None)
    parser.add_argument("--min-country-regions", type=int, default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    model_cfg = config["analysis"]["socioeconomic_model"]
    alpha = (
        float(args.alpha)
        if args.alpha is not None
        else float(model_cfg["alpha"])
    )
    min_regions = (
        int(args.min_country_regions)
        if args.min_country_regions is not None
        else int(model_cfg["min_country_regions"])
    )

    data = prepare_data(config, min_regions)
    processed = configured_path(config, "processed")

    print(
        f"Rows used: {len(data):,}; "
        f"countries: {data['country'].nunique()}; "
        f"repeated NUTS3 groups: {data['NUTS_ID'].nunique():,}"
    )
    print(
        "Maximal model: log_pop * log_gdp_pc * year_centered * C(country) "
        "+ (1 | NUTS_ID)"
    )
    print(f"Deletion threshold alpha = {alpha:g}")

    terms = all_hierarchical_terms()
    formula = formula_from_terms(terms)
    current, optimizer, fit_warnings = fit_model(formula, data)

    history = []
    deletion_step = 0

    # Crawley-style simplification: highest order to lowest order.
    for order in range(len(TOKENS), 0, -1):
        while True:
            candidates = removable_terms(terms, order)
            if not candidates:
                break

            tests = []
            for candidate in candidates:
                reduced_terms = set(terms)
                reduced_terms.remove(candidate)
                reduced_formula = formula_from_terms(reduced_terms)
                reduced, reduced_optimizer, reduced_warnings = fit_model(
                    reduced_formula,
                    data,
                )
                lr, df, p = lrt(current, reduced)
                tests.append(
                    {
                        "candidate": candidate,
                        "reduced_terms": reduced_terms,
                        "reduced_formula": reduced_formula,
                        "reduced_result": reduced,
                        "optimizer": reduced_optimizer,
                        "warnings": " | ".join(reduced_warnings),
                        "lr": lr,
                        "df": df,
                        "p": p,
                    }
                )

            # Delete the least significant removable term if non-significant.
            chosen = max(
                tests,
                key=lambda item: (
                    -np.inf if np.isnan(item["p"]) else item["p"]
                ),
            )
            delete = (
                np.isfinite(chosen["p"])
                and chosen["p"] >= alpha
            )

            for test in tests:
                history.append(
                    {
                        "deletion_step": deletion_step,
                        "interaction_order": order,
                        "current_formula": formula_from_terms(terms),
                        "term_tested": term_string(test["candidate"]),
                        "lr_statistic": test["lr"],
                        "df_removed": test["df"],
                        "p_value": test["p"],
                        "selected_for_deletion": (
                            delete
                            and test["candidate"] == chosen["candidate"]
                        ),
                        "decision": (
                            "delete"
                            if delete
                            and test["candidate"] == chosen["candidate"]
                            else "retain_at_this_step"
                        ),
                        "reduced_formula": test["reduced_formula"],
                        "optimizer": test["optimizer"],
                        "warnings": test["warnings"],
                    }
                )

            if not delete:
                print(
                    f"Order {order}: all {len(candidates)} removable term(s) "
                    f"significant at alpha={alpha:g}; move to lower order."
                )
                break

            deletion_step += 1
            removed_name = term_string(chosen["candidate"])
            print(
                f"Delete {removed_name}: LR={chosen['lr']:.4f}, "
                f"df={chosen['df']}, p={chosen['p']:.6g}"
            )
            terms = chosen["reduced_terms"]
            formula = chosen["reduced_formula"]
            current = chosen["reduced_result"]
            optimizer = chosen["optimizer"]
            fit_warnings = chosen["warnings"]

    final_formula = formula_from_terms(terms)
    marginal_r2, conditional_r2 = nakagawa_r2(current)
    random_var = float(np.asarray(current.cov_re)[0, 0])
    residual_var = float(current.scale)
    icc = random_var / (random_var + residual_var)

    fixed = pd.DataFrame(
        {
            "term": current.fe_params.index,
            "estimate": current.fe_params.values,
            "std_error": current.bse_fe.values,
        }
    )
    fixed["z"] = fixed["estimate"] / fixed["std_error"]
    fixed["p_value"] = 2 * norm.sf(np.abs(fixed["z"]))

    retained = pd.DataFrame(
        [
            {
                "term": term_string(term),
                "order": len(term),
                "required_by_hierarchy": any(
                    term < other for other in terms
                ),
            }
            for term in sorted(terms, key=term_key)
        ]
    )

    history_df = pd.DataFrame(history)
    history_path = processed / "socioeconomic_model_simplification.csv"
    fixed_path = processed / "socioeconomic_minimal_model_fixed_effects.csv"
    retained_path = processed / "socioeconomic_minimal_model_terms.csv"
    summary_path = processed / "socioeconomic_minimal_model_summary.txt"

    history_df.to_csv(history_path, index=False)
    fixed.to_csv(fixed_path, index=False)
    retained.to_csv(retained_path, index=False)

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("CRAWLEY-STYLE MINIMUM ADEQUATE MIXED MODEL\n")
        handle.write("=" * 72 + "\n\n")
        handle.write(f"Alpha: {alpha}\n")
        handle.write(f"Observations: {len(data)}\n")
        handle.write(f"Countries: {data['country'].nunique()}\n")
        handle.write(f"NUTS3 groups: {data['NUTS_ID'].nunique()}\n")
        handle.write(f"Optimizer: {optimizer}\n")
        handle.write(f"Formula: {final_formula}\n\n")
        handle.write(current.summary().as_text())
        handle.write("\n\n")
        handle.write(f"Nakagawa marginal R2: {marginal_r2:.8f}\n")
        handle.write(f"Nakagawa conditional R2: {conditional_r2:.8f}\n")
        handle.write(f"NUTS3 random-intercept variance: {random_var:.8f}\n")
        handle.write(f"Residual variance: {residual_var:.8f}\n")
        handle.write(f"NUTS3 ICC: {icc:.8f}\n")

    print("\nMINIMUM ADEQUATE MODEL")
    print("=" * 90)
    print(final_formula)
    print()
    print(
        f"logLik={current.llf:.6f}; AIC={current.aic:.3f}; "
        f"BIC={current.bic:.3f}"
    )
    print(
        f"Marginal R2={marginal_r2:.6f}; "
        f"conditional R2={conditional_r2:.6f}; "
        f"NUTS3 ICC={icc:.6f}"
    )
    print()
    print("Retained terms:")
    print(retained.to_string(index=False))
    print()
    print(f"Wrote: {history_path}")
    print(f"Wrote: {fixed_path}")
    print(f"Wrote: {retained_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
