#!/usr/bin/env python3
"""Descriptive spatial dilution diagnostics for five documented lighting cases.

Uses existing locally extracted EOG V1 VCMSLCFG 75N060W monthly avg_rade9h
and cf_cvg TIFFs; NO further downloads or annual/PRC recalculation.
Strata are defined ONLY by *pre-intervention reference-year* radiance or a
fixed geographic bounding-box quadrant, never by the observed event change.
Treatments are NOT inferred from these strata or from declining pixels.

    python src/local_viirs_spatial_dilution.py --self-test
    python src/local_viirs_spatial_dilution.py --cases portsmouth bristol
    python src/local_viirs_spatial_dilution.py
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd

from local_viirs_five_case_analysis import (
    CASES, SECONDARY, build_inventory, grid_month, load_boundary,
    ROOT, UNIT, _outline,
)

OUT = ROOT / "data/processed/local_lighting_interventions/spatial_dilution"
FIG = ROOT / "figures/local_lighting_interventions/spatial_dilution"
BRIGHTNESS_LABELS = (
    "lowest_50pct", "50_to_80pct", "80_to_95pct", "top_5pct"
)
GEO_LABELS = ("NW", "NE", "SW", "SE")
# Months in the outcome year preceding any independently documented rollout.
# None for Thanet: Dec 2016-May 2017 was a staged transition.
PRE_MONTHS = {
    "echirolles": (4, 5, 6),
    "portsmouth": (1, 2, 3, 4, 5),
    "coventry": (1, 2, 3, 4),
    "bristol": (1, 2),
    "thanet": (),
}


def strict_pair(grids, outcome_year, reference_years, month, min_cf=1):
    """One same-calendar-month, strictly paired native-pixel delta.

    To be comparable each cell needs valid radiance and >=min_cf observations
    in the outcome AND EACH reference year for that same calendar month.
    """
    current = grids[f"{outcome_year}-{month:02d}"]
    rr = [grids[f"{y}-{month:02d}"] for y in reference_years]
    a = current["rad"]
    bstack = np.stack([z["rad"] for z in rr])
    coverage = np.stack([current["cf"]] + [z["cf"] for z in rr])
    common = (
        np.isfinite(a) & np.all(np.isfinite(bstack), axis=0)
        & np.all(np.isfinite(coverage) & (coverage >= min_cf), axis=0)
    )
    b = np.mean(bstack, axis=0)
    baseline = np.full(a.shape, np.nan, dtype=np.float32)
    current_valid = np.full(a.shape, np.nan, dtype=np.float32)
    baseline[common] = b[common]
    current_valid[common] = a[common]
    return baseline, current_valid, current_valid - baseline


def period_maps(grids, year, references, months, *, min_cf=1, support=.75):
    """Mean month-matched baseline/outcome/Δ on a fixed within-period grid.

    An eligible cell has valid outcome-and-all-reference data in at least
    ceil(support * number of specified calendar months). The same observed
    monthly cell support is used for baseline, outcome, and their difference.
    """
    pairs = [strict_pair(grids, year, references, m, min_cf=min_cf)
             for m in months]
    base = np.stack([p[0] for p in pairs])
    outcome = np.stack([p[1] for p in pairs])
    n = np.sum(np.isfinite(base), axis=0)
    target = max(1, math.ceil(len(months) * support))
    def mean(z):
        dst = np.full(z.shape[1:], np.nan, dtype=np.float32)
        np.divide(np.nansum(z, axis=0), n, out=dst, where=n >= target)
        return dst
    b, a = mean(base), mean(outcome)
    return b, a, a - b, n


def baseline_strata(grids, reference_years, months, inside):
    """Fixed brightness ranks based solely on reference-year VIIRS inputs.

    Uses the mean reference radiance across ALL pre+post calendar months
    prespecified for the *primary* event, independently of outcome year.
    Requires at least half of those reference observations available.
    Brightness groups are rank-defined: 50%, 30%, 15%, highest 5%.
    Tie handling is deterministic; groups are descriptive, not exposures.
    """
    stacks = np.stack([
        grids[f"{y}-{m:02d}"]["rad"] for y in reference_years for m in months
    ])
    count = np.isfinite(stacks).sum(axis=0)
    valid = inside & (count >= math.ceil(stacks.shape[0] / 2))
    avg = np.full(stacks.shape[1:], np.nan, dtype=np.float32)
    np.divide(np.nansum(stacks, axis=0), count, out=avg, where=valid)
    flat = np.flatnonzero(valid)
    if len(flat) < 12:
        raise ValueError("Fewer than 12 baseline cells; insufficient for strata")
    order = flat[np.argsort(avg.flat[flat], kind="stable")]
    ranks = np.full(avg.size, -1, dtype=np.int8)
    for i, segment in enumerate(np.array_split(np.arange(len(order)), len(order))):
        pass  # reserve no indices; explicit ranked thresholds below
    n = len(order)
    for i, index in enumerate(order):
        q = (i + .5) / n
        ranks[index] = 0 if q < .50 else 1 if q < .80 else 2 if q < .95 else 3
    return avg, ranks.reshape(avg.shape)


def geographic_quadrants(grids, months, bounds):
    """Fixed NW/NE/SW/SE city bounding-box partition in EPSG:4326.

    No VIIRS outcome can change this assignment. Source monthly EOG gridding
    is EPSG:4326; refuse to silently mislabel another projection.
    """
    grid = next(iter(grids.values()))
    if grid["crs"].to_epsg() != 4326:
        raise ValueError("Geographic quadrants require EPSG:4326 EOG raster")
    transform = grid["transform"]
    yy, xx = np.indices(grid["rad"].shape)
    lon = transform.c + transform.a * (xx + .5) + transform.b * (yy + .5)
    lat = transform.f + transform.d * (xx + .5) + transform.e * (yy + .5)
    minx, miny, maxx, maxy = bounds
    west = lon < (minx + maxx) / 2
    north = lat >= (miny + maxy) / 2
    out = np.full(xx.shape, -1, dtype=np.int8)
    inside = grid["inside"]
    out[inside & west & north] = 0
    out[inside & ~west & north] = 1
    out[inside & west & ~north] = 2
    out[inside & ~west & ~north] = 3
    return out


def group_summary(key, design, bin_type, labels, bins, base, after, diff,
                  baseline_level, pre_delta, min_cf):
    """Report exact arithmetic spatial decomposition of the municipal mean.

    Group contribution = sum paired-pixel ΔR in group / total paired pixels.
    These contributions sum to the matched-pixel mean ΔR, not to the original
    variable-valid-citymean estimator. Negative/positive shares are kept
    separately because net cancellation may be substantial.
    """
    eligible = np.isfinite(diff) & (bins >= 0)
    n_all = int(eligible.sum())
    if not n_all:
        raise ValueError(f"{key} {design}: no paired pixels")
    bsum = float(np.nansum(base[eligible]))
    losses = float(np.maximum(-diff[eligible], 0).sum())
    gains = float(np.maximum(diff[eligible], 0).sum())
    rows = []
    memberships = [("all_pixels", eligible)]
    memberships += [(lab, eligible & (bins == i))
                    for i, lab in enumerate(labels)]
    for label, mask in memberships:
        n = int(mask.sum())
        if not n:
            continue
        b, a, d = base[mask], after[mask], diff[mask]
        matched_pre = pre_delta[mask] if pre_delta is not None else np.array([])
        valid_pre = np.isfinite(matched_pre)
        post_on_common = d[valid_pre] if pre_delta is not None else np.array([])
        pre_on_common = matched_pre[valid_pre] if pre_delta is not None else np.array([])
        local_loss = float(np.maximum(-d, 0).sum())
        local_gain = float(np.maximum(d, 0).sum())
        mean_b = float(b.mean())
        median_b_anchor = float(np.nanmedian(baseline_level[mask]))
        rows.append({
            "case": key, "contrast": design, "minimum_cf_per_month": min_cf,
            "stratification": bin_type, "group": label,
            "n_paired_pixels": n, "n_city_paired_pixels": n_all,
            "fraction_city_pixels": n / n_all,
            "reference_mean": mean_b, "outcome_mean": float(a.mean()),
            "mean_delta": float(d.mean()), "median_delta": float(np.median(d)),
            "percent_change": (100 * float(d.mean()) / mean_b
                               if mean_b > .05 else np.nan),
            "fraction_pixels_negative": float((d < 0).mean()),
            "baseline_brightness_group_median": median_b_anchor,
            "share_baseline_radiance": (
                float(b.sum()) / bsum if bsum > 0 else np.nan),
            "weighted_net_delta_to_city_mean": float(d.sum() / n_all),
            "negative_delta_magnitude": local_loss,
            "positive_delta_magnitude": local_gain,
            "share_all_decline_magnitude": (
                local_loss / losses if losses > 0 else np.nan),
            "share_all_increase_magnitude": (
                local_gain / gains if gains > 0 else np.nan),
            "n_pre_and_post_paired_pixels": int(valid_pre.sum()),
            "pre_mean_delta_common": (
                float(pre_on_common.mean()) if valid_pre.any() else np.nan),
            "post_mean_delta_common": (
                float(post_on_common.mean()) if valid_pre.any() else np.nan),
            "post_minus_pre_common": (
                float((post_on_common - pre_on_common).mean())
                if valid_pre.any() else np.nan),
            "caveat": ("Baseline-defined groups describe concentration only; "
                       "no streetlight treatment polygon has been established."),
        })
    return rows


def make_map(ax, arr, ref, title, *, cmap, norm=None, vmin=None, vmax=None):
    l, b, r, t = ref["extent"]
    image = ax.imshow(
        arr, extent=(l, r, b, t), origin="upper", interpolation="nearest",
        cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
    _outline(ax, ref["geometry"])
    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=8)
    return image


def plot_diagnostics(key, name, geom, strata, sectors, baseline, pre, post,
                     counts, design):
    stem = f"{key}_{design}"
    FIG.mkdir(parents=True, exist_ok=True)
    ref = baseline["grid"]
    anchor = baseline["anchor"]
    pre_delta = pre["delta"] if pre is not None else None
    post_delta = post["delta"]
    diff = (post_delta - pre_delta if pre_delta is not None else
            np.full(post_delta.shape, np.nan))
    vals = np.r_[post_delta[np.isfinite(post_delta)],
                 diff[np.isfinite(diff)]]
    scale = max(float(np.nanquantile(np.abs(vals), .98)), 1e-3)
    norm = TwoSlopeNorm(-scale, vcenter=0, vmax=scale)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10),
                             constrained_layout=True)
    palette = ListedColormap(["#d4e6fa", "#92b7d8", "#477ca4", "#183b69"])
    cmapnorm = BoundaryNorm(np.arange(-.5, 4.5), palette.N)
    x = np.where(ref["inside"] & (strata >= 0), strata, np.nan)
    m = make_map(axes[0, 0], x, ref, "Pre-event reference brightness rank",
                 cmap=palette, norm=cmapnorm)
    cb = fig.colorbar(m, ax=axes[0, 0], shrink=.78, ticks=np.arange(4))
    cb.ax.set_yticklabels(["bottom 50%", "50–80%", "80–95%", "top 5%"])
    make_map(axes[0, 1], pre_delta if pre_delta is not None
             else np.full(post_delta.shape, np.nan), ref,
             "Pre-event seasonal anomaly" if pre is not None else
             "No clean within-year pre-event period", norm=norm,
             cmap="coolwarm")
    m = make_map(axes[1, 0], post_delta, ref, "Post-event seasonal anomaly",
                 norm=norm, cmap="coolwarm")
    make_map(axes[1, 1], diff, ref,
             "Post minus pre anomaly (same pixels)"
             if pre is not None else "Not estimable without a clean pre-period",
             norm=norm, cmap="coolwarm")
    fig.colorbar(m, ax=[axes[0, 1], axes[1, 0], axes[1, 1]],
                 shrink=.75, label=f"Δ radiance ({UNIT})")
    fig.suptitle(
        f"{name}: spatial concentration, fixed pre-event strata\n"
        "Strictly paired same-month VIIRS pixels; not treatment assignment")
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / f"{stem}_spatial_anomalies{ext}", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5),
                             constrained_layout=True)
    valid_post = np.isfinite(post_delta)
    sector_map = np.where(valid_post, sectors, np.nan)
    m1 = make_map(axes[0], sector_map, ref, "Fixed geographic quadrants",
                  cmap=palette, norm=cmapnorm)
    cb = fig.colorbar(m1, ax=axes[0], ticks=np.arange(4))
    cb.ax.set_yticklabels(GEO_LABELS)
    m2 = make_map(axes[1], np.where(ref["inside"], counts, np.nan), ref,
                  "Post-period matched monthly observations",
                  cmap="cividis", vmin=0, vmax=post["nmonths"])
    fig.colorbar(m2, ax=axes[1],
                 label=f"Supported months (out of {post['nmonths']})")
    fig.suptitle(f"{name}: fixed geography and valid-month support")
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / f"{stem}_geographic_support{ext}", dpi=190)
    plt.close(fig)


def plot_strata(key, design, rows):
    df = pd.DataFrame(rows)
    df = df[(df.stratification == "baseline_rank") &
            (df.group != "all_pixels") &
            (df.minimum_cf_per_month == 1)]
    if df.empty:
        return
    order = list(BRIGHTNESS_LABELS)
    df = df.set_index("group").reindex(order).dropna(subset=["mean_delta"])
    fig, axes = plt.subplots(2, 1, figsize=(10, 8),
                             constrained_layout=True)
    x = np.arange(len(df))
    axes[0].bar(x, df["mean_delta"])
    axes[0].axhline(0, color="black", linewidth=.8)
    axes[0].set_ylabel(f"Mean paired-cell Δ ({UNIT})")
    axes[1].bar(x, df["weighted_net_delta_to_city_mean"])
    axes[1].axhline(0, color="black", linewidth=.8)
    axes[1].set_ylabel(f"Contribution to municipal mean Δ ({UNIT})")
    axes[1].set_xticks(x, df.index.tolist())
    for ax in axes:
        ax.grid(axis="y", alpha=.3)
    fig.suptitle(
        f"{key.title()}: baseline brightness bands — {design}\n"
        "Bands defined from reference years only; NOT a treatment map")
    for ext in (".png", ".pdf"):
        fig.savefig(FIG / f"{key}_{design}_brightness_decomposition{ext}",
                    dpi=190)
    plt.close(fig)


def run_case(key, index):
    case = CASES[key]
    geom, geo_label, geo_source = load_boundary(key, case)
    primary_months = set(case.compare_months) | set(PRE_MONTHS[key])
    specs = [("primary", case.treatment_year,
              case.reference_years, case.compare_months)]
    if key in SECONDARY:
        y, refs, months, _ = SECONDARY[key]
        specs.append(("secondary", y, refs, months))
    required = set()
    for _, y, refs, months in specs:
        for m in set(months) | set(PRE_MONTHS[key]):
            required.add(f"{y}-{m:02d}")
            required.update(f"{r}-{m:02d}" for r in refs)
    for y in case.reference_years:
        required.update(f"{y}-{m:02d}" for m in primary_months)

    print(f"\n{key.upper()} | geography: {geo_label} | "
          f"needed months: {len(required)}", flush=True)
    grids, sig = {}, None
    for ym in sorted(required):
        grid = grid_month(geom, *index[ym])
        signature = (grid["crs"], grid["transform"], grid["rad"].shape)
        if sig is None:
            sig = signature
        elif sig != signature:
            raise ValueError(f"{key} {ym}: native grid does not align")
        grids[ym] = grid
    refgrid = next(iter(grids.values()))
    anchor, strata = baseline_strata(
        grids, case.reference_years, sorted(primary_months),
        refgrid["inside"])
    sectors = geographic_quadrants(grids, primary_months, geom.bounds)
    rows = []
    for tag, year, refs, months in specs:
        design = f"{tag}_{year}"
        # Strict same-cell monthly comparison at cf_cvg>=1,3,5, all references.
        for threshold in (1, 3, 5):
            b, a, d, n = period_maps(
                grids, year, refs, months, min_cf=threshold)
            pre = None
            if PRE_MONTHS[key]:
                pb, pa, pd, pn = period_maps(
                    grids, year, refs, PRE_MONTHS[key], min_cf=threshold)
                pre = {"delta": pd, "nmonth": pn}
            for family, labels, groups in (
                ("baseline_rank", BRIGHTNESS_LABELS, strata),
                ("geographic_quadrant", GEO_LABELS, sectors),
            ):
                these = group_summary(
                    key, design, family, labels, groups, b, a, d,
                    anchor, pre["delta"] if pre is not None else None,
                    threshold)
                rows.extend(these)
            if threshold == 1:
                plot_diagnostics(
                    key, case.name, geom, strata, sectors,
                    {"grid": refgrid, "anchor": anchor}, pre,
                    {"delta": d, "nmonths": len(months)}, n,
                    design)
                plot_strata(key, design, rows)
                baseline_city = float(np.mean(b[np.isfinite(d)]))
                mean_city = float(np.mean(d[np.isfinite(d)]))
                print(f"  {design}: strict matched ΔR={mean_city:+.3f}; "
                      f"reference mean={baseline_city:.3f}, "
                      f"support={int(np.isfinite(d).sum())} cells", flush=True)
    return rows, {"case": key, "boundary": geo_source,
                  "baseline_reference_years": ",".join(
                      map(str, case.reference_years)),
                  "baseline_strata_calendar_months": ",".join(
                      map(str, sorted(primary_months))),
                  "baseline_cells": int((strata >= 0).sum())}


def self_test():
    small = {}
    for y in (2019, 2020, 2021):
        for m in (1, 2):
            a = np.array([[1, 2, 3, 4],
                          [5, 6, 7, 8],
                          [9, 10, 11, 12]], dtype=np.float32)
            if y == 2021:
                a = a - np.array([[1, 0, 1, 0],
                                  [0, 2, 0, 2],
                                  [0, 3, 0, 3]], dtype=np.float32)
            small[f"{y}-{m:02d}"] = {
                "rad": a, "cf": np.full(a.shape, 8, dtype=np.float32),
                "inside": np.ones(a.shape, dtype=bool)}
    anchor, strata = baseline_strata(
        small, (2019, 2020), (1, 2),
        small["2019-01"]["inside"])
    assert int((strata >= 0).sum()) == 12
    assert len(set(strata.ravel())) == 4
    b, a, delta, n = period_maps(
        small, 2021, (2019, 2020), (1, 2))
    assert np.isclose(delta[1, 1], -2)
    assert int(n[1, 1]) == 2
    rows = group_summary(
        "synthetic", "test", "baseline_rank",
        BRIGHTNESS_LABELS, strata, b, a, delta, anchor, None, 1)
    group = [r for r in rows if r["group"] != "all_pixels"]
    allrow = rows[0]
    assert np.isclose(sum(r["weighted_net_delta_to_city_mean"]
                          for r in group), allrow["mean_delta"])
    assert np.isclose(sum(r["share_baseline_radiance"]
                          for r in group), 1)
    _, _, d, n = period_maps(
        small, 2021, (2019, 2020), (1, 2), min_cf=9)
    assert not np.isfinite(d).any()
    print("Self-test passed: fixed baseline ranks, strict common cells, "
          "cf_cvg masks, additive contributions.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=list(CASES),
                        default=list(CASES))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    index, missing = build_inventory()
    if missing:
        raise RuntimeError(
            f"Missing {len(missing)} EOG VCMSLCFG monthly input layers; "
            "run src/inventory_monthly_viirs_acquisitions.py")
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    allrows, metadata = [], []
    for key in args.cases:
        these, meta = run_case(key, index)
        allrows.extend(these)
        metadata.append(meta)
    results = pd.DataFrame(allrows)
    results.to_csv(OUT / "spatial_stratified_contrasts.csv", index=False)
    pd.DataFrame(metadata).to_csv(
        OUT / "baseline_spatial_strata_metadata.csv", index=False)
    subset = results[
        (results.minimum_cf_per_month == 1) &
        (results.stratification == "baseline_rank")]
    subset.to_csv(OUT / "brightness_strata_summary.csv", index=False)
    print(f"\nDone: {len(args.cases)} cities; {len(results)} stratified rows.")
    print(f"CSV: {OUT}\nMaps and bar charts: {FIG}")
    print("No intervention group inferred from the outcome; "
          "predefined brightness/geographic groups only. "
          "No annual VIIRS/PRC files modified.")


if __name__ == "__main__":
    main()
