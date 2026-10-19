#!/usr/bin/env python3
"""Refit R SSasymp on provisional site-matched PRC and observation-year EOG VIIRS.

Requires: numpy, pandas, scipy, matplotlib. No remote downloads.
Inputs: digitized original plotted 18 PRC observations (approximate radiances,
provisional site IDs) and the actual locally extracted site annual EOG VIIRS CSV.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def predict(x, z):
    yf, y0, log_alpha = z
    return yf + (y0 - yf) * np.exp(-np.exp(log_alpha) * np.asarray(x, float))


def fit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 4 or not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("Fit needs at least four finite pairs")
    candidates = [
        least_squares(lambda z: predict(x, z) - y,
                      [yf, max(y), np.log(k)], max_nfev=4000)
        for yf in (0.10, 0.15, 0.20)
        for k in (0.04, 0.10, 0.25)
    ]
    z = min(candidates, key=lambda a: np.dot(a.fun, a.fun)).x
    residual = y - predict(x, z)
    return {
        "z": z, "y_f": z[0], "y_0": z[1],
        "log_alpha": z[2], "k": np.exp(z[2]), "n": len(x),
        "R2": 1 - sum(residual**2) / sum((y - y.mean())**2),
        "RMSE": float(np.sqrt(np.mean(residual**2))),
    }


def heldout_site_rmse(d, xcol):
    errors = []
    for site_id in d.site_id_candidate.unique():
        test = d.site_id_candidate.eq(site_id)
        trained = fit(d.loc[~test, xcol], d.loc[~test, "prc_sigg"])
        errors.extend((d.loc[test, "prc_sigg"].to_numpy() -
                       predict(d.loc[test, xcol], trained["z"])).tolist())
    return float(np.sqrt(np.mean(np.square(errors))))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignment", type=Path,
                    default=Path("data/reference/prc_alignment_from_pdf_provisional.csv"))
    ap.add_argument("--site-comparison", type=Path,
                    default=Path("data/processed/prc_site_radiance_2018_2019/prc_site_annual_radiance_comparison.csv"))
    ap.add_argument("--output-dir", type=Path,
                    default=Path("data/processed/prc_fit_observation_year"))
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(a.alignment)
    sites = pd.read_csv(a.site_comparison)
    assert len(obs) == 18 and obs.sigg_index.is_unique
    assert len(sites) == 14 and sites.site_id.is_unique
    s = sites[["site_id", "field_date", "eog_2018_pixel_corrected",
               "eog_2019_pixel_corrected"]].copy()
    d = obs.merge(s, left_on="site_id_candidate", right_on="site_id",
                  how="left", validate="many_to_one")
    years = pd.to_datetime(d.field_date).dt.year
    d["radiance_corrected_field_year"] = np.where(
        years.eq(2018), d.eog_2018_pixel_corrected,
        np.where(years.eq(2019), d.eog_2019_pixel_corrected, np.nan))
    matched = d[d.radiance_corrected_field_year.notna()].copy()
    assert len(matched) == 16 and matched.site_id_candidate.nunique() == 14
    unique = matched[matched.site_mapping_status.eq("unique")].copy()
    assert len(unique) == 12

    original = "original_plot_radiance_approx"
    updated = "radiance_corrected_field_year"
    specs = [
        ("Original 18, approximate", d, original),
        ("Original matched 16", matched, original),
        ("Corrected EOG matched 16", matched, updated),
        ("Original unambiguous 12", unique, original),
        ("Corrected EOG unambiguous 12", unique, updated),
    ]
    results = {}
    for name, frame, xcol in specs:
        info = fit(frame[xcol], frame.prc_sigg)
        info["site_holdout_RMSE"] = (
            heldout_site_rmse(frame, xcol)
            if frame.site_id_candidate.notna().all() else np.nan)
        results[name] = info
    pd.DataFrame([
        {"model": name, **{key: value for key, value in info.items() if key != "z"}}
        for name, info in results.items()
    ]).to_csv(a.output_dir / "prc_fit_model_comparison.csv", index=False)

    for key, xcol, name in [
        ("original18", original, "Original 18, approximate"),
        ("original16", original, "Original matched 16"),
        ("corrected16", updated, "Corrected EOG matched 16"),
    ]:
        d["predicted_" + key] = predict(d[xcol], results[name]["z"])
    d.to_csv(a.output_dir / "prc_fit_input_predictions.csv", index=False)

    grid = np.linspace(0, 65, 500)
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2),
                            sharex=True, sharey=True, layout="constrained")
    for ax, name, xcol, title in [
        (axs[0], "Original matched 16", original, "A   Original plotted radiance"),
        (axs[1], "Corrected EOG matched 16", updated,
         "B   Corrected EOG in observation year"),
    ]:
        m = results[name]
        for yr, sym in [(2018, "s"), (2019, "o")]:
            rows = matched[pd.to_datetime(matched.field_date).dt.year.eq(yr)]
            ax.scatter(rows[xcol], rows.prc_sigg, s=48, marker=sym,
                       label=f"{yr} field observations", zorder=4)
        ax.plot(grid, predict(grid, m["z"]), lw=2, label="SSasymp fit")
        ax.text(.98, .96, f"n=16; sites=14\nR²={m['R2']:.3f}\nRMSE={m['RMSE']:.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=10)
        ax.set_title(title, loc="left")
        ax.set_xlim(-1, 65)
        ax.set_ylim(.09, .51)
        ax.set_xlabel("Radiance (nW cm$^{-2}$ sr$^{-1}$)")
        ax.grid(alpha=.2)
        ax.legend(loc="upper right", bbox_to_anchor=(1, .79), frameon=False, fontsize=9)
    axs[0].set_ylabel("Polarization-receptor contrast (PRC)")
    fig.suptitle("Original asymptotic PRC function: radiance predictor comparison")
    for extension in ("png", "pdf"):
        fig.savefig(a.output_dir / ("prc_fit_annual_radiance_comparison." + extension),
                    dpi=240, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for name in ["Original 18, approximate", "Original matched 16",
                 "Corrected EOG matched 16"]:
        ax.plot(grid, predict(grid, results[name]["z"]), label=name)
    ax.set(xlim=(0, 65), ylim=(.1, .49), xlabel="Radiance (nW cm$^{-2}$ sr$^{-1}$)",
           ylabel="Predicted PRC", title="Comparison of fitted transfer functions")
    ax.grid(alpha=.2)
    ax.legend()
    for extension in ("png", "pdf"):
        fig.savefig(a.output_dir / ("prc_fit_transfer_curves." + extension),
                    dpi=240, bbox_inches="tight")
    plt.close(fig)

    lines = [
        "Unweighted R SSasymp equivalent: PRC(R)=yf+(y0-yf)*exp(-exp(log_alpha)*R)",
        "Site matching provisional; 2/18 observations lack identifiable locations.",
        "Same 16 observations used for main predictor comparison.",
        "CV grouped by proposed site ID, keeping apparent replicates together.",
    ]
    for name, m in results.items():
        lines.append(
            f"{name}: n={m['n']} y0={m['y_0']:.6f} yf={m['y_f']:.6f} "
            f"k={m['k']:.6f} R2={m['R2']:.6f} RMSE={m['RMSE']:.6f} "
            f"site-held-out RMSE={m['site_holdout_RMSE']:.6f}")
    (a.output_dir / "prc_fit_summary.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print("Saved:", a.output_dir)


if __name__ == "__main__":
    main()
