#!/usr/bin/env python3
"""Post-extraction diagnostic checks for the *observed* five-case EOG VIIRS results.

This analyses only outputs from local_viirs_five_case_analysis.py; it neither
downloads files nor reads the raw monthly GeoTIFFs. No fitting, p-values, PRC
recalculation or causal attribution is performed. Outputs in separate
five_case_robustness subdirectories.

Run:
    python src/local_viirs_five_case_robustness.py --self-test
    python src/local_viirs_five_case_robustness.py
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/processed/local_lighting_interventions/five_case_monthly"
OUT = ROOT / "data/processed/local_lighting_interventions/five_case_robustness"
FIG = ROOT / "figures/local_lighting_interventions/five_case_robustness"
UNIT = "nW cm$^{-2}$ sr$^{-1}$"


@dataclass(frozen=True)
class Spec:
    outcome: int
    references: tuple[int, ...]
    after: tuple[int, ...]
    before: tuple[int, ...]
    transition: tuple[int, ...]
    placebos: tuple[tuple[int, tuple[int, ...]], ...]
    caveat: str


SPECS = {
    "echirolles": Spec(2022, (2020, 2021), (10, 11, 12),
                      (4, 5, 6), (7, 8, 9),
                      ((2021, (2019, 2020)), (2020, (2018, 2019))),
                      "Before=Apr-Jun to avoid Jul-Sep earlier switch-offs; 2023 LED changes."),
    "portsmouth": Spec(2020, (2018, 2019), (9, 10, 11, 12),
                       (1, 2, 3, 4, 5), (6, 7, 8),
                       ((2019, (2017, 2018)), (2018, (2016, 2017))),
                       "Earlier 2018/19 LED rollout contaminates reference/placebos."),
    "coventry": Spec(2024, (2022, 2023), (8, 9, 10, 11, 12),
                     (1, 2, 3, 4), (5, 6, 7),
                     ((2023, (2021, 2022)), (2022, (2020, 2021))),
                     "Jan-Apr 2024 predates May-Jul switch-off; Nov/Dec may have low cf_cvg."),
    "bristol": Spec(2023, (2020, 2021, 2022), tuple(range(4, 13)),
                    (1, 2), (3,),
                    ((2022, (2019, 2020, 2021)),),
                    "Only Jan-Feb 2023 fully precede staged rollout; post is partly treated."),
    "thanet": Spec(2017, (2015, 2016), (6, 7, 8, 9, 10, 11),
                   (), (1, 2, 3, 4, 5),
                   ((2016, (2015,)),),
                   "2017 Jan-May belongs to staged installation, so no within-2017 pre."),
}

VARS = (
    ("all_observed", "mean_radiance"),
    ("cf_cvg_3plus", "mean_cf3"),
    ("cf_cvg_5plus", "mean_cf5"),
)


def check_input(key):
    path = SOURCE / f"{key}_2015_2024_monthly.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run src/local_viirs_five_case_analysis.py first.")
    df = pd.read_csv(path)
    needed = ("year", "calendar_month", "mean_radiance", "median_cf_cvg",
              "valid_percent", "mean_cf3", "mean_cf5")
    absent = [c for c in needed if c not in df]
    if absent:
        raise ValueError(f"{key}: missing required columns: {absent}")
    if df.duplicated(["year", "calendar_month"]).any():
        raise ValueError(f"{key}: duplicate year/calendar-month records")
    if len(df) != 120:
        raise ValueError(f"{key}: expected 120 observations, got {len(df)}")
    return df.set_index(["year", "calendar_month"]).sort_index()


def compare(panel, year, reference_years, months, var="mean_radiance",
            min_median_cf=None):
    """Compare same-calendar-month municipal means on eligible month support.

    Missing or low-coverage months are excluded *as whole monthly units*, not
    treated as zero. Data values are city means, not fixed matched pixels.
    """
    details = []
    for m in months:
        current = panel.loc[(year, m)]
        references = panel.loc[[(y, m) for y in reference_years]]
        obs = float(current[var])
        base_values = references[var].astype(float)
        reference = float(base_values.mean()) if base_values.notna().all() else np.nan
        obs_cf = float(current["median_cf_cvg"])
        ref_cf = float(references["median_cf_cvg"].min())
        include = (
            np.isfinite(obs) and np.isfinite(reference) and
            (min_median_cf is None or
             (obs_cf >= min_median_cf and ref_cf >= min_median_cf))
        )
        details.append({
            "month": m, "observed": obs, "reference": reference,
            "delta": obs-reference if include else np.nan,
            "included": include, "outcome_cf": obs_cf,
            "min_reference_cf": ref_cf,
            "outcome_valid_percent": float(current["valid_percent"]),
        })
    d = pd.DataFrame(details)
    used = d[d.included]
    if used.empty:
        return {"delta_abs": np.nan, "delta_pct": np.nan, "n_months": 0,
                "observed_mean": np.nan, "reference_mean": np.nan,
                "median_cf_outcome": np.nan, "valid_outcome_mean": np.nan}, d
    observed = used["observed"].mean()
    reference = used["reference"].mean()
    return {
        "delta_abs": float(observed-reference),
        "delta_pct": float(100*(observed-reference)/reference)
                     if reference > .05 else np.nan,
        "n_months": len(used),
        "observed_mean": float(observed),
        "reference_mean": float(reference),
        "median_cf_outcome": float(used.outcome_cf.median()),
        "valid_outcome_mean": float(used.outcome_valid_percent.mean()),
    }, d


def one_case(key, panel, spec):
    rows = []
    base_results = {}
    for label, col in VARS:
        for refs_label, refs in [
            ("all_reference_years", spec.references),
            *[(f"reference_{year}", (year,)) for year in spec.references],
        ]:
            for filter_name, threshold in (
                ("all_months", None), ("median_cf5_months_only", 5),
                ("median_cf10_months_only", 10)):
                result, _ = compare(panel, spec.outcome, refs, spec.after,
                                    col, min_median_cf=threshold)
                rows.append({
                    "case": key, "scenario": "primary", "year": spec.outcome,
                    "metric": label, "reference_type": refs_label,
                    "reference_years": ",".join(map(str, refs)),
                    "month_selection": filter_name,
                    "months_requested": ",".join(map(str, spec.after)),
                    "omitted_month": "", **result,
                })
                if (label == "all_observed" and
                    refs_label == "all_reference_years" and
                    filter_name == "all_months"):
                    base_results = result
    # Leave-one-month-out: is a large effect driven by a single month?
    for drop in spec.after:
        months = tuple(m for m in spec.after if m != drop)
        result, _ = compare(panel, spec.outcome, spec.references, months)
        rows.append({
            "case": key, "scenario": "leave_one_month_out",
            "year": spec.outcome, "metric": "all_observed",
            "reference_type": "all_reference_years",
            "reference_years": ",".join(map(str, spec.references)),
            "month_selection": "all_except_one",
            "months_requested": ",".join(map(str, months)),
            "omitted_month": drop, **result,
        })

    pre_results = None
    if spec.before:
        pre_results, _ = compare(panel, spec.outcome,
                                 spec.references, spec.before)
        prepost = {
            "case": key, "year": spec.outcome,
            "reference_years": ",".join(map(str, spec.references)),
            "pre_months": ",".join(map(str, spec.before)),
            "post_months": ",".join(map(str, spec.after)),
            "pre_anomaly": pre_results["delta_abs"],
            "post_anomaly": base_results["delta_abs"],
            "post_minus_pre_anomaly": (
                base_results["delta_abs"] - pre_results["delta_abs"]),
            "pre_months_used": pre_results["n_months"],
            "post_months_used": base_results["n_months"],
            "note": spec.caveat,
        }
    else:
        prepost = {
            "case": key, "year": spec.outcome,
            "reference_years": ",".join(map(str, spec.references)),
            "pre_months": "", "post_months": ",".join(map(str, spec.after)),
            "pre_anomaly": np.nan, "post_anomaly": base_results["delta_abs"],
            "post_minus_pre_anomaly": np.nan, "pre_months_used": 0,
            "post_months_used": base_results["n_months"],
            "note": spec.caveat + " No within-year pre-period comparator.",
        }

    placebos = []
    for year, refs in spec.placebos:
        result, _ = compare(panel, year, refs, spec.after)
        placebos.append({
            "case": key, "placebo_year": year,
            "reference_years": ",".join(map(str, refs)),
            "months": ",".join(map(str, spec.after)), **result,
            "note": spec.caveat + " Placebo is an earlier calendar comparison, "
                    "not necessarily a genuinely untreated control.",
        })
    return rows, prepost, placebos


def plot_case(key, panel, spec):
    # Use same baseline each month rather than a seasonally confounded
    # raw pre/post difference.
    months = np.arange(1, 13)
    obs = np.array([panel.loc[(spec.outcome, m), "mean_radiance"] for m in months])
    baseline = np.array([
        np.mean([panel.loc[(y, m), "mean_radiance"] for y in spec.references])
        for m in months])
    obs_cf = np.array([panel.loc[(spec.outcome, m), "median_cf_cvg"] for m in months])
    ref_cf = np.array([
        np.mean([panel.loc[(y, m), "median_cf_cvg"] for y in spec.references])
        for m in months])
    fig, axes = plt.subplots(3, 1, figsize=(11, 9),
                             constrained_layout=True, sharex=True)
    axes[0].plot(months, obs, "-o", label=str(spec.outcome), markersize=4)
    axes[0].plot(months, baseline, "-o",
                 label="Mean reference ("+",".join(map(str, spec.references))+")",
                 markersize=4)
    axes[0].set_ylabel(f"Mean monthly VIIRS ({UNIT})")
    axes[0].legend()
    axes[1].plot(months, obs-baseline, "-o", markersize=4)
    axes[1].axhline(0, color="black", linewidth=.8)
    axes[1].set_ylabel(f"Same-month anomaly ({UNIT})")
    axes[2].plot(months, obs_cf, "-o", label=str(spec.outcome))
    axes[2].plot(months, ref_cf, "-o", label="Reference")
    axes[2].set_ylabel("Median cloud-free observations")
    axes[2].legend()
    for ax in axes:
        if spec.transition:
            ax.axvspan(min(spec.transition)-.45,
                       max(spec.transition)+.45, color="gray",
                       alpha=.18, label="Transition")
        ax.grid(alpha=.2)
    axes[-1].set(xlim=(.6, 12.4), xticks=months, xlabel="Calendar month")
    fig.suptitle(
        f"{key.title()}: {spec.outcome} vs same-calendar-month references. "
        "Shaded: transition.\nNo causal effect inferred.")
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / f"{key}_event_anomalies{ext}", dpi=190)
    plt.close(fig)


def self_test():
    rows = []
    for y in (2019, 2020, 2021):
        for m in range(1, 13):
            rows.append({
                "year": y, "calendar_month": m,
                "mean_radiance": (10. if y < 2021 else 8.),
                "mean_cf3": (10. if y < 2021 else 8.),
                "mean_cf5": (10. if y < 2021 else 8.),
                "median_cf_cvg": 9, "valid_percent": 100,
            })
    panel = pd.DataFrame(rows).set_index(["year", "calendar_month"])
    result, _ = compare(panel, 2021, (2019, 2020), (9, 10, 11, 12))
    assert result["n_months"] == 4
    assert np.isclose(result["delta_abs"], -2)
    assert np.isclose(result["delta_pct"], -20)
    result2, _ = compare(panel, 2021, (2019, 2020), (9, 10),
                         min_median_cf=10)
    assert result2["n_months"] == 0
    print("Self-test passed: same-calendar-month contrasts and coverage rule.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--cases", choices=list(SPECS), nargs="+",
                   default=list(SPECS))
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rows, preposts, placebos = [], [], []
    for key in args.cases:
        panel = check_input(key)
        r, p, pl = one_case(key, panel, SPECS[key])
        rows.extend(r)
        preposts.append(p)
        placebos.extend(pl)
        plot_case(key, panel, SPECS[key])
        print(f"{key}: pre-event Δ={p['pre_anomaly']:+.2f}, "
              f"post-event Δ={p['post_anomaly']:+.2f}"
              if np.isfinite(p["pre_anomaly"])
              else f"{key}: no defensible within-outcome-year pre-period",
              flush=True)
    pd.DataFrame(rows).to_csv(OUT / "reference_and_cloudfree_sensitivity.csv",
                              index=False)
    pd.DataFrame(preposts).to_csv(OUT / "pre_post_anomaly_differentials.csv",
                                  index=False)
    pd.DataFrame(placebos).to_csv(OUT / "earlier_year_placebo_comparisons.csv",
                                  index=False)
    print(f"Results: {OUT}\nFigures: {FIG}")
    print("Descriptive checks ONLY. Cf thresholding changes pixel support; "
          "placebo earlier years may contain LED changes; no causal model fitted.")


if __name__ == "__main__":
    main()
