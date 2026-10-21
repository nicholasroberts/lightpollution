#!/usr/bin/env python3
"""Test whether historical ISS colour changes covary with corrected VIIRS dimming.

Read-only join of existing paired ISS NUTS2016 RGB comparisons and fixed NUTS2016
2014/2020/2024 annual radiance. No refitting, source overwrites, PRC correction,
or assertion of annual LED adoption timing. Results are DESCRIPTIVE; the ISS
comparison is 2012–13 versus an irregular 2014–20 composite and therefore
cannot directly measure a specific 2014–24 spectral transition.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from project_config import PROJECT_ROOT as ROOT

ISS = ROOT / "data/processed/iss_nuts3_spectral_shift.csv"
VIIRS = ROOT / "data/processed/fixed_nuts2016_viirs_2014_2024.parquet"
OUT = ROOT / "data/processed/iss_spectral_vs_corrected_viirs_2014_2024"
FIG = ROOT / "figures/iss_spectral_vs_corrected_viirs_2014_2024"
YEARS = (2014, 2020, 2024)
WINDOWS = ((2014, 2020), (2014, 2024), (2020, 2024))
COLOURS = ("delta_gr_median", "delta_bg_median")
LOG15 = float(np.log10(1.5))


def spearman(x, y):
    a, b = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 4 or len(np.unique(a)) < 2 or len(np.unique(b)) < 2:
        return float("nan")
    return float(spearmanr(a, b).statistic)


def prepare(iss, rad, min_pixels):
    required_iss = {"NUTS_ID", "CNTR_CODE", "n_common_valid_pixels",
                    "delta_gr_median", "delta_bg_median"}
    required_rad = {"NUTS_ID", "year", "radiance_mean_corrected"}
    if required_iss - set(iss) or required_rad - set(rad):
        raise ValueError("Missing source columns: ISS "
                         f"{required_iss-set(iss)}; VIIRS {required_rad-set(rad)}")
    iss = iss.copy()
    iss.NUTS_ID = iss.NUTS_ID.astype(str)
    if iss.NUTS_ID.duplicated().any():
        raise ValueError("Duplicate ISS NUTS3 identifier")
    rad = rad[rad.year.isin(YEARS)].copy()
    rad.NUTS_ID = rad.NUTS_ID.astype(str)
    if rad.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate VIIRS year/region")
    wide = rad.pivot(index="NUTS_ID", columns="year",
                     values="radiance_mean_corrected")
    if not set(YEARS).issubset(set(wide.columns)):
        raise ValueError("Missing 2014/2020/2024 VIIRS years")
    wide = wide[list(YEARS)].rename(columns=lambda year: f"radiance_{year}")
    d = iss.merge(wide.reset_index(), on="NUTS_ID", how="inner",
                  validate="one_to_one")
    d["n_common_valid_pixels"] = pd.to_numeric(
        d.n_common_valid_pixels, errors="coerce")
    for col in [*COLOURS, *[f"radiance_{year}" for year in YEARS]]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["eligible_iss_and_radiance"] = (
        d.n_common_valid_pixels.ge(min_pixels)
        & np.isfinite(d[[*COLOURS, *[f"radiance_{y}" for y in YEARS]]]
                      .to_numpy(dtype=float)).all(axis=1)
        & d[[f"radiance_{y}" for y in YEARS]].gt(0).all(axis=1)
    )
    for first, last in WINDOWS:
        d[f"delta_log10_radiance_{first}_{last}"] = (
            np.log10(d[f"radiance_{last}"].where(d[f"radiance_{last}"] > 0))
            - np.log10(d[f"radiance_{first}"].where(d[f"radiance_{first}"] > 0)))
        d[f"delta_radiance_{first}_{last}"] = (
            d[f"radiance_{last}"] - d[f"radiance_{first}"])
    # Published empirical G_ISS / VIIRS proportional to 1.5^(G_ISS / R_ISS)
    # is a bounded equal-GREEN-OUTPUT thought experiment. Do NOT present as
    # a recovered real flux, actual VIIRS correction or LED conversion rate.
    d["illustrative_equal_green_spectral_delta_log10_viirs"] = (
        -LOG15 * d.delta_gr_median)
    d["illustrative_2014_2020_non_spectral_residual"] = (
        d.delta_log10_radiance_2014_2020
        - d.illustrative_equal_green_spectral_delta_log10_viirs)
    return d


def associations(d):
    results = []
    eligible = d[d.eligible_iss_and_radiance]
    for n_min in (20, 50):
        q = eligible[eligible.n_common_valid_pixels.ge(n_min)]
        for country, part in [("ALL", q), *list(q.groupby("CNTR_CODE"))]:
            for first, last in WINDOWS:
                outcome = f"delta_log10_radiance_{first}_{last}"
                for spectral in COLOURS:
                    results.append({
                        "scope": "pooled_regions" if country == "ALL" else "within_country",
                        "country": country, "min_common_iss_pixels": n_min,
                        "n": len(part), "radiance_window": f"{first}-{last}",
                        "spectral_metric": spectral,
                        "spearman_rho": spearman(part[spectral], part[outcome]),
                        "median_spectral_change": (
                            float(part[spectral].median()) if len(part) else np.nan),
                        "median_viirs_log_change": (
                            float(part[outcome].median()) if len(part) else np.nan),
                        "n_viirs_down": int(part[outcome].lt(0).sum()),
                        "n_viirs_down_and_spectral_up": int(
                            (part[outcome].lt(0) & part[spectral].gt(0)).sum()),
                    })
    return pd.DataFrame(results)


def national(d):
    """Equal-country contrast of regional medians (NOT pooled NUTS3 regression)."""
    rows = []
    for country, z in d[d.eligible_iss_and_radiance].groupby("CNTR_CODE"):
        rec = {"country": country, "n": len(z)}
        for s in COLOURS:
            rec[s] = float(z[s].median())
        for f, t in WINDOWS:
            rec[f"median_delta_log10_radiance_{f}_{t}"] = float(
                z[f"delta_log10_radiance_{f}_{t}"].median())
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("country")


def country_associations(c):
    rows = []
    for nmin in (5, 10, 20):
        q = c[c.n >= nmin]
        for f, t in WINDOWS:
            target = f"median_delta_log10_radiance_{f}_{t}"
            for spectral in COLOURS:
                rows.append({
                    "min_regions_per_country": nmin, "n_countries": len(q),
                    "radiance_window": f"{f}-{t}",
                    "spectral_metric": spectral,
                    "spearman_rho_equal_country_weight": spearman(q[spectral], q[target])
                })
    return pd.DataFrame(rows)


def sign_table(d):
    """Signs only, with uncertainty around near-zero spectral change explicit."""
    rows = []
    q = d[d.eligible_iss_and_radiance]
    for cutoff in (0.0, 0.02):
        for first, last in WINDOWS:
            y = q[f"delta_log10_radiance_{first}_{last}"]
            for spectral in COLOURS:
                pos = q[spectral].gt(cutoff)
                neg = q[spectral].lt(-cutoff)
                neutral = ~(pos | neg)
                for class_name, mask in [("bluer_or_whiter", pos),
                                         ("less_blue_or_white", neg),
                                         ("near_zero_or_uncertain", neutral)]:
                    z = y[mask]
                    rows.append({"radiance_window": f"{first}-{last}",
                                 "spectral_metric": spectral, "delta_ratio_cutoff": cutoff,
                                 "spectral_group": class_name, "n": len(z),
                                 "n_dimmed_viirs": int(z.lt(0).sum()),
                                 "pct_dimmed_viirs": 100 * float(z.lt(0).mean())
                                     if len(z) else np.nan,
                                 "median_log10_viirs_change": float(z.median())
                                     if len(z) else np.nan})
    return pd.DataFrame(rows)


def figures(d, countries, path):
    path.mkdir(parents=True, exist_ok=True)
    d = d[d.eligible_iss_and_radiance]
    fig, ax = plt.subplots(2, 2, figsize=(13, 11), layout="constrained")
    for j, spectral in enumerate(COLOURS):
        for i, (start, end) in enumerate(((2014, 2020), (2014, 2024))):
            a = ax[j, i]
            response = f"delta_log10_radiance_{start}_{end}"
            a.axhline(0, lw=.8, color="gray")
            a.axvline(0, lw=.8, color="gray")
            a.scatter(d[spectral], d[response], s=8, alpha=.32)
            rho = spearman(d[spectral], d[response])
            a.set(title=f"{start}–{end} VIIRS; pooled regional rho={rho:.2f}, n={len(d)}",
                  xlabel=f"ISS paired colour-ratio change: {spectral}",
                  ylabel="log10(corrected annual VIIRS endpoint ratio)")
    fig.suptitle("Pre-2013 vs 2014–20 ISS spectral shift against 2014–20 / 2014–24 VIIRS\n"
                 "Epochs are not contemporaneously aligned; correlation ≠ spectral attribution")
    fig.savefig(path/"iss_spectral_shift_vs_corrected_viirs_scatter.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for i, spectral in enumerate(COLOURS):
        q = countries[countries.n >= 10]
        outcome = "median_delta_log10_radiance_2014_2020"
        ax[i].axvline(0, lw=.8, color="gray")
        ax[i].axhline(0, lw=.8, color="gray")
        ax[i].scatter(q[spectral], q[outcome], s=27)
        for row in q.itertuples(index=False):
            if row.country in ("FR", "DE", "IT", "NL", "IE", "UK", "BE", "RO"):
                ax[i].annotate(row.country, (getattr(row, spectral), getattr(row, outcome)),
                               xytext=(3, 3), textcoords="offset points", fontsize=8)
        rho = spearman(q[spectral], q[outcome])
        ax[i].set(title=f"Country medians, n={len(q)}, rho={rho:.2f}",
                  xlabel=spectral, ylabel="Median NUTS3 Δlog10 VIIRS 2014–2020")
    fig.suptitle("One country per observation; >=10 spectral-valid regions")
    fig.savefig(path/"iss_spectral_shift_vs_viirs_country_medians.pdf", bbox_inches="tight")
    plt.close(fig)


def self_test():
    s = pd.DataFrame({"NUTS_ID": ["DE001", "DE002", "FR001", "FR002"],
                      "CNTR_CODE": ["DE", "DE", "FR", "FR"],
                      "delta_gr_median": [0., .1, -.1, .1],
                      "delta_bg_median": [.1, .1, 0., -.1],
                      "n_common_valid_pixels": [100, 30, 50, 10]})
    r = pd.DataFrame([{"NUTS_ID": nid, "year": year,
                       "radiance_mean_corrected": val * ((1.+delta)**((year-2014)/6.))}
                      for nid, val, delta in
                      [("DE001", 2., -.1), ("DE002", 3., -.2),
                       ("FR001", 1., .05), ("FR002", 5., -.1)]
                      for year in YEARS])
    d = prepare(s, r, 20)
    assert len(d) == 4 and int(d.eligible_iss_and_radiance.sum()) == 3
    assert np.isclose(float(d.loc[d.NUTS_ID.eq("DE001"),
                                   "delta_log10_radiance_2014_2020"].iloc[0]),
                      np.log10(.9))
    assert np.isclose(float(d.loc[d.NUTS_ID.eq("DE002"),
                                   "illustrative_equal_green_spectral_delta_log10_viirs"].iloc[0]),
                      -.1 * LOG15)
    assert len(associations(d)) > 0 and len(sign_table(d)) > 0
    print("ISS/VIIRS paired-endpoint synthetic QA: PASS")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--iss", type=Path, default=ISS)
    p.add_argument("--viirs", type=Path, default=VIIRS)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--figures", type=Path, default=FIG)
    p.add_argument("--min-iss-pixels", type=int, default=20)
    a = p.parse_args()
    if a.min_iss_pixels < 1:
        p.error("--min-iss-pixels must be >=1")
    self_test()
    if a.self_test:
        return
    if not a.iss.exists():
        raise FileNotFoundError(f"Missing {a.iss}; run existing ISS spectral validation first")
    if not a.viirs.exists():
        raise FileNotFoundError(f"Missing {a.viirs}; use locked zero-offset VIIRS annual product")
    iss = pd.read_csv(a.iss, dtype={"NUTS_ID": str, "CNTR_CODE": str})
    rad = pd.read_parquet(a.viirs)
    d = prepare(iss, rad, a.min_iss_pixels)
    if d.eligible_iss_and_radiance.sum() < 50:
        raise ValueError("Too few valid ISS/VIIRS NUTS3 pairs; inspect source coverage")
    within = associations(d)
    c = national(d)
    between = country_associations(c)
    signs = sign_table(d)
    a.out.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out/"iss_viirs_matched_nuts3.csv", index=False)
    within.to_csv(a.out/"iss_viirs_within_country_associations.csv", index=False)
    c.to_csv(a.out/"iss_viirs_country_medians.csv", index=False)
    between.to_csv(a.out/"iss_viirs_between_country_associations.csv", index=False)
    signs.to_csv(a.out/"iss_viirs_sign_sensitivity.csv", index=False)
    figures(d, c, a.figures)
    (a.out/"iss_viirs_analysis_manifest.json").write_text(json.dumps({
        "iss": str(a.iss), "viirs": str(a.viirs), "eligible_regions":
        int(d.eligible_iss_and_radiance.sum()),
        "iss_ratio_epochs": ["2012-2013", "2014-2020 irregular mosaic"],
        "annual_viirs_endpoints": list(YEARS),
        "minimum_iss_common_pixels": a.min_iss_pixels,
        "2021_snow": "No 2021 used; therefore results do not depend on whether a NUTS3 2021 snow-adjusted radiance table has been published.",
        "pixel_support": "ISS paired pre/post common footprint, different from VIIRS NUTS3 whole-polygon mean",
        "illustrative_equal_green": "-log10(1.5)*paired_delta_G_over_R; not a calibrated light-emission correction",
        "inference": "Descriptive associations only; no spatially-independent p-values and no year-resolved LED timing",
        "source_table_status": "READ ONLY; original ISS/VIIRS/PRC sources not rewritten",
    }, indent=2) + "\n")
    print(f"Joined {int(d.eligible_iss_and_radiance.sum())} ISS/VIIRS NUTS3 regions")
    print(within[(within.country.isin(["ALL", "FR", "UK", "IT", "NL", "DE", "IE", "BE", "RO"]))
                 & within.radiance_window.eq("2014-2020")
                 & within.spectral_metric.eq("delta_gr_median")
                 & within.min_common_iss_pixels.eq(20)]
          [["country", "n", "spearman_rho", "median_spectral_change",
            "median_viirs_log_change"]].to_string(index=False))
    print("Output:", a.out)


if __name__ == "__main__":
    main()
