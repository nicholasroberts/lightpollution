#!/usr/bin/env python3
"""Publication Figure 4A–D: regional demographic vs PRC sign contrasts.

READ ONLY: uses frozen 1,273 matched NUTS2016 Jan-1 density/VIIRS/PRC
endpoints for 2014–2020; final locked spider/moth pixel-derived PRC
annual tables for THREE illustrative trajectories; and existing Jan-1
regional annual density. No PRC fitting, source correction, new
population reconstruction or population modelling.

A. Spider and moth 2x2 sign-classification matrices (counts, fractions).
B. Two NUTS3 maps of same categorical classification (spider and moth).
C. Regional category proportions within six previously chosen countries.
D. Three regional 2014–2020 change trajectories: FR, TR, NL333.
   FR/TR chosen deterministically as central (rather than extreme)
   examples of specified sign combinations; user may supply actual
   fixed NUTS3 IDs with --fr-id and --tr-id before final paper selection.

Default sign classification is strictly positive vs strictly negative;
exact-zero changes are neutral (never coded as decline). Optional neutral
deadbands are recorded in the output and are display/category sensitivity
choices, not reanalysis. Outputs editable SVG+PDF and audited CSV/JSON.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT, load_config
from figure_style import apply_figure_style, figure_size, save_figure, SERIES_COLOURS

ROOT = PROJECT_ROOT
ENDPOINTS = ROOT / "data/processed/jan1_population_radiance_prc_country_nl333/matched_region_change_2014_2020.csv"
ANNUAL_POP = ROOT / "data/processed/nuts2016_jan1_population_density_sensitivity/regional_annual_density_2014_2020_full_audit.csv"
ANNUAL_PRC = {
    taxon: ROOT / "data/processed" / f"historical_prc_2014_2024_{taxon}_snow_corrected" / "prc_nuts3_annual_cumulative.csv"
    for taxon in ("spider", "moth")
}
OUTPUT = ROOT / "figures/publication/figure4_regional_contrasts"
YEARS = tuple(range(2014, 2021))
MATCHED_N = 1273
POP = "change_density_jan1_fixed_land_2014_2020"
PRC = {
    "spider": "cumulative_delta_spider_prc_2020",
    "moth": "cumulative_delta_moth_prc_2020",
}
COUNTRIES = ("TR", "DE", "FR", "IE", "IT", "NL")
COUNTRY_NAMES = {
    "TR": "Türkiye", "DE": "Germany", "FR": "France",
    "IE": "Ireland", "IT": "Italy", "NL": "Netherlands",
    "UK": "United Kingdom",
}
# Order retained in all charts, map keys, output tables and manifests.
CATEGORIES = (
    "growth_prc_loss", "growth_prc_gain",
    "decline_prc_loss", "decline_prc_gain", "neutral"
)
COLORS = {
    "growth_prc_loss": "#D55E00",
    "growth_prc_gain": "#009E73",
    "decline_prc_loss": "#7B61A8",
    "decline_prc_gain": "#0072B2",
    "neutral": "#989898",
}
LABELS = {
    "growth_prc_loss": "Population ↑ / PRC ↓",
    "growth_prc_gain": "Population ↑ / PRC ↑",
    "decline_prc_loss": "Population ↓ / PRC ↓",
    "decline_prc_gain": "Population ↓ / PRC ↑",
    "neutral": "Neutral (zero or within threshold)",
}
EXTENT = (-12.0, 34.0, 46.0, 73.0)  # Includes eastern Türkiye.
LABEL_POP = "Change in Jan-1 density from 2014 (people km⁻²)"
LABEL_PRC = "Change in PRC from 2014 (dimensionless)"


def resolve(path: Path | str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def frozen_endpoints(path: Path, *, expected_n: int = MATCHED_N) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing frozen 2014–2020 Jan1/PRC source {path}. "
            "Do not substitute the old 1,317-unit official-density series."
        )
    d = pd.read_csv(path, dtype={"NUTS_ID": str, "country": str})
    required = {
        "NUTS_ID", "country", POP, *PRC.values(),
        "density_jan1_fixed_land_2014", "density_jan1_fixed_land_2020",
        "fixed_land_area_km2", "radiance_2014", "radiance_2020",
    }
    if missing := required-set(d):
        raise ValueError(f"Frozen joined endpoints missing {sorted(missing)}")
    if len(d) != expected_n or not d.NUTS_ID.is_unique or d.NUTS_ID.isna().any():
        raise ValueError(f"Expected {expected_n} distinct NUTS2016 IDs; got {len(d)}")
    if not d.country.str.fullmatch(r"[A-Z]{2}").all():
        raise ValueError("Invalid country labels in frozen endpoints")
    for col in required-{"NUTS_ID", "country"}:
        d[col] = pd.to_numeric(d[col], errors="coerce")
        if not np.isfinite(d[col]).all():
            raise ValueError(f"Non-finite frozen {col}")
    if (d.fixed_land_area_km2 <= 0).any():
        raise ValueError("Invalid fixed land area")
    if (d[["density_jan1_fixed_land_2014", "density_jan1_fixed_land_2020"]] <= 0).any().any():
        raise ValueError("Invalid Jan-1 density endpoints")
    if not np.allclose(
        d[POP],
        d["density_jan1_fixed_land_2020"]-d["density_jan1_fixed_land_2014"],
        rtol=1e-10, atol=1e-8
    ):
        raise ValueError("Frozen population change is not the endpoint difference")
    if "NL333" not in set(d.NUTS_ID) or d.loc[d.NUTS_ID.eq("NL333"), "country"].iloc[0] != "NL":
        raise ValueError("Documented NL333 case missing from frozen cohort")
    return d.sort_values(["country", "NUTS_ID"]).reset_index(drop=True)


def categories(pop: np.ndarray, prc: np.ndarray,
               pop_deadband: float, prc_deadband: float) -> np.ndarray:
    """Every observation assigned exactly once; no fabricated sign at zero."""
    p = np.asarray(pop, dtype=float)
    y = np.asarray(prc, dtype=float)
    if not np.isfinite(p).all() or not np.isfinite(y).all():
        raise ValueError("Classification requires finite signed changes")
    if pop_deadband < 0 or prc_deadband < 0:
        raise ValueError("Neutral thresholds must be >= 0")
    result = np.full(p.shape, "neutral", dtype=object)
    up, down = p > pop_deadband, p < -pop_deadband
    gain, loss = y > prc_deadband, y < -prc_deadband
    result[up & loss] = "growth_prc_loss"
    result[up & gain] = "growth_prc_gain"
    result[down & loss] = "decline_prc_loss"
    result[down & gain] = "decline_prc_gain"
    return result


def classified(d: pd.DataFrame, pop_deadband: float, prc_deadband: float):
    x = d.copy()
    for taxon in PRC:
        x[f"category_{taxon}"] = categories(
            x[POP].to_numpy(float), x[PRC[taxon]].to_numpy(float),
            pop_deadband, prc_deadband
        )
        assert x[f"category_{taxon}"].isin(CATEGORIES).all()
        assert len(x[f"category_{taxon}"]) == len(x)
    return x


def count_table(d: pd.DataFrame) -> pd.DataFrame:
    groups = []
    for taxon in PRC:
        for cat in CATEGORIES:
            n = int(d[f"category_{taxon}"].eq(cat).sum())
            groups.append({
                "taxon": taxon, "category": cat, "n_nuts3": n,
                "percent_matched_nuts3": 100*n/len(d)
            })
    counts = pd.DataFrame(groups)
    if not counts.groupby("taxon").n_nuts3.sum().eq(len(d)).all():
        raise ValueError("Sign classifications did not partition cohort")
    return counts


def plot_4a(counts: pd.DataFrame, out: Path, n: int):
    """Paired four-quadrant counts with neutral cases reported beneath."""
    fig, axs = plt.subplots(1, 2, figsize=(figure_size("slide")[0], 4.25))
    table = (
        ("decline_prc_gain", "growth_prc_gain"),
        ("decline_prc_loss", "growth_prc_loss"),
    )
    for ax, taxon in zip(axs, PRC):
        data = counts[counts.taxon.eq(taxon)].set_index("category")
        for i, row in enumerate(table):
            for j, cat in enumerate(row):
                patch = Rectangle((j, i), 1, 1, facecolor=COLORS[cat],
                                  edgecolor="white", linewidth=2)
                ax.add_patch(patch)
                count = int(data.loc[cat, "n_nuts3"])
                pct = float(data.loc[cat, "percent_matched_nuts3"])
                ax.text(j+.5, i+.5, f"{count:,} regions\n{pct:.1f}% of cohort",
                        ha="center", va="center", color="white",
                        weight="bold", fontsize=9)
        neutral = int(data.loc["neutral", "n_nuts3"])
        ax.set_xlim(0, 2)
        ax.set_ylim(2, 0)
        ax.set_aspect("equal")
        ax.set_xticks([.5, 1.5], ["Population decreased", "Population increased"])
        ax.set_yticks([.5, 1.5], ["PRC increased", "PRC decreased"])
        ax.tick_params(axis="both", length=0, pad=8)
        for label in ax.get_xticklabels():
            label.set_fontsize(8)
        ax.set_title(f"{taxon.capitalize()} (n={n:,})", loc="left", fontsize=10)
        ax.text(.02, 1.015, f"Neutral: {neutral} regions", transform=ax.transAxes,
                ha="left", fontsize=8, color="#555555")
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("4A  |  Signs of regional population and PRC change, 2014–2020",
                 x=.03, ha="left", fontsize=11)
    fig.subplots_adjust(left=.15, right=.99, bottom=.15, top=.82, wspace=.3)
    save_figure(fig, out)
    plt.close(fig)


def load_display_geometry(config: dict, tolerance_m: float):
    if not np.isfinite(tolerance_m) or tolerance_m < 0:
        raise ValueError("Display simplification must be nonnegative and finite")
    nuts = load_nuts(config)
    if len(nuts) != 1522 or nuts.NUTS_ID.duplicated().any():
        raise ValueError("The original fixed NUTS2016 1,522 polygons are required")
    envelope = box(*EXTENT)
    clipped = gpd.clip(nuts, envelope)
    clipped = clipped.loc[~clipped.geometry.is_empty].to_crs("EPSG:3035")
    if tolerance_m:
        clipped["geometry"] = clipped.geometry.simplify(
            tolerance=tolerance_m, preserve_topology=True
        )
    projected_extent = gpd.GeoSeries([envelope], crs="EPSG:4326").to_crs("EPSG:3035")
    return (clipped[["NUTS_ID", "CNTR_CODE", "geometry"]],
            tuple(projected_extent.total_bounds))


def plot_4b(geom: gpd.GeoDataFrame, bounds: tuple,
            d: pd.DataFrame, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 6.6))
    fig.subplots_adjust(left=.02, right=.99, bottom=.17, top=.90, wspace=.045)
    country_outline = geom.dissolve(by="CNTR_CODE").boundary
    data = d[["NUTS_ID", "category_spider", "category_moth"]]
    g = geom.merge(data, on="NUTS_ID", how="left", validate="one_to_one")
    x0, y0, x1, y1 = bounds
    mapped = {}
    for ax, taxon in zip(axes, PRC):
        key = f"category_{taxon}"
        colors = g[key].map(COLORS).fillna("#E8E8E8")
        g.plot(ax=ax, color=colors, edgecolor="none", linewidth=0)
        country_outline.plot(ax=ax, color="#484848", linewidth=.28)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(f"{taxon.capitalize()} PRC", loc="left", fontsize=10)
        mapped[taxon] = int(g[key].notna().sum())
    fig.suptitle("4B  |  Regional sign combinations across Europe, 2014–2020",
                 x=.025, ha="left", fontsize=11)
    legend = [
        Patch(facecolor=COLORS[cat], label=LABELS[cat])
        for cat in CATEGORIES
    ] + [Patch(facecolor="#E8E8E8", label="Outside matched population cohort")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               fontsize=8, bbox_to_anchor=(.5, .012))
    save_figure(fig, out, map_figure=True)
    plt.close(fig)
    return mapped


def plot_4c(d: pd.DataFrame, countries: tuple[str, ...], out: Path):
    """Country category fractions retain the *same* country sample in each taxon."""
    fig, ax = plt.subplots(figsize=(figure_size("slide")[0], 6.25))
    country_counts = []
    ylabels, yvalues = [], []
    for code in countries:
        subset = d[d.country.eq(code)]
        if not len(subset):
            raise ValueError(f"Country absent in frozen 2014–2020 source: {code}")
        for taxon in ("spider", "moth"):
            ylabels.append(f"{COUNTRY_NAMES.get(code, code)} · {taxon.capitalize()}")
            yvalues.append((code, taxon, subset))
    yy = np.arange(len(yvalues))
    for yi, (code, taxon, subset) in enumerate(yvalues):
        x0 = 0.
        for cat in CATEGORIES:
            n = int(subset[f"category_{taxon}"].eq(cat).sum())
            fraction = 100*n/len(subset)
            ax.barh(yi, fraction, left=x0, height=.72, color=COLORS[cat],
                    edgecolor="white", linewidth=.28)
            x0 += fraction
            country_counts.append({
                "country": code, "taxon": taxon, "category": cat,
                "n_regions": n, "country_total_nuts3": len(subset),
                "percentage_of_country_matched_nuts3": fraction
            })
        if not np.isclose(x0, 100., atol=1e-10):
            raise ValueError("Country category proportions did not sum to 100")
        ax.text(101, yi, f"n={len(subset)}", fontsize=8, va="center")
    ax.set_yticks(yy, ylabels)
    ax.invert_yaxis()
    ax.set_xlim(0, 114)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Share of matched NUTS3 regions (%)")
    ax.set_title("4C  |  Country differences in regional sign categories, 2014–2020",
                 loc="left", fontsize=10)
    ax.grid(axis="x", alpha=.13)
    ax.set_axisbelow(True)
    fig.legend(
        handles=[Patch(facecolor=COLORS[cat], label=LABELS[cat]) for cat in CATEGORIES],
        ncol=2, frameon=False, loc="lower center", fontsize=8,
        bbox_to_anchor=(.53, -.005)
    )
    fig.subplots_adjust(left=.27, right=.99, top=.92, bottom=.20)
    save_figure(fig, out)
    plt.close(fig)
    return pd.DataFrame(country_counts)


def central_example(d: pd.DataFrame, code: str, prc_sign: int,
                    *, pop_deadband: float, prc_deadband: float) -> tuple[str, int]:
    """Choose *central* illustrative region, not most extreme PRC outlier.

    FR: Δdensity positive, spider ΔPRC positive.
    TR: Δdensity positive, spider ΔPRC negative.
    Both selections are deterministic conditional summaries, NOT assertions
    of representativeness or mechanistic causation.
    """
    if prc_sign not in (-1, 1):
        raise ValueError("Expected PRC sign -1/+1")
    sub = d[d.country.eq(code)]
    mask = sub[POP] > pop_deadband
    mask &= (sub[PRC["spider"]]*prc_sign) > prc_deadband
    cand = sub.loc[mask].copy()
    if len(cand) == 0:
        raise ValueError(
            f"No {code} region matches selection; give a valid --{code.lower()}-id "
            "or adjust explicitly labelled neutral sensitivity thresholds."
        )
    vpop = cand[POP].to_numpy(float)
    vprc = cand[PRC["spider"]].to_numpy(float)
    scale_pop = max(np.quantile(vpop, .75)-np.quantile(vpop, .25), 1e-8)
    scale_prc = max(np.quantile(vprc, .75)-np.quantile(vprc, .25), 1e-8)
    cand["_distance"] = (
        np.abs(vpop-np.median(vpop))/scale_pop +
        np.abs(vprc-np.median(vprc))/scale_prc
    )
    cand = cand.sort_values(["_distance", "NUTS_ID"])
    return str(cand.iloc[0].NUTS_ID), len(cand)


def example_ids(d: pd.DataFrame, args):
    out = {}
    choices = {"FR": (args.fr_id, 1), "TR": (args.tr_id, -1),
               "NL": (args.nl_id, None)}
    selection = {}
    for code, (user_id, prc_sign) in choices.items():
        if user_id:
            selected, rule, eligible = str(user_id), "explicit NUTS3 ID", None
        elif code == "NL":
            selected, rule, eligible = "NL333", "previously investigated NL333", None
        else:
            selected, eligible = central_example(
                d, code, prc_sign,
                pop_deadband=args.pop_deadband, prc_deadband=args.prc_deadband
            )
            rule = "minimum IQR-scaled distance to within-sign-group medians"
        matches = d[d.NUTS_ID.eq(selected)]
        if len(matches) != 1 or str(matches.iloc[0].country) != code:
            raise ValueError(f"Selected {selected} is not a unique matched {code} NUTS3")
        out[code] = selected
        selection[code] = {"NUTS_ID": selected, "rule": rule,
                           "eligible_regions_for_selection": eligible,
                           "signed_population_change": float(matches.iloc[0][POP]),
                           "signed_spider_prc_change": float(matches.iloc[0][PRC["spider"]]),
                           "signed_moth_prc_change": float(matches.iloc[0][PRC["moth"]])}
    return out, selection


def trajectories(ids: dict[str, str], endpoints: pd.DataFrame,
                 pop_source: Path, prc_sources: dict[str, Path]) -> pd.DataFrame:
    """Read the original full annual observations for three selected regions."""
    if not pop_source.is_file():
        raise FileNotFoundError(f"Missing frozen annual Jan1 audit: {pop_source}")
    selected = set(ids.values())
    pop = pd.read_csv(pop_source, dtype={"NUTS_ID": str})
    required = {"NUTS_ID", "year", "density_jan1_fixed_land",
                "reconstruction_eligible_year"}
    if missing := required-set(pop):
        raise ValueError(f"Annual Jan1 audit missing columns {sorted(missing)}")
    pop = pop[pop.NUTS_ID.isin(selected)&pop.year.isin(YEARS)].copy()
    if len(pop) != 3*len(YEARS) or pop.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Selected annual Jan1 density series incomplete")
    eligibility = pop.reconstruction_eligible_year.astype(str).str.lower()
    if not eligibility.isin({"true", "1"}).all():
        raise ValueError("One of the frozen annual demographic records is ineligible")
    pop["density_jan1_fixed_land"] = pd.to_numeric(
        pop.density_jan1_fixed_land, errors="coerce"
    )
    if not np.isfinite(pop.density_jan1_fixed_land).all() or (
        pop.density_jan1_fixed_land <= 0
    ).any():
        raise ValueError("Invalid annual Jan1 density values")
    result = pop[["NUTS_ID", "year", "density_jan1_fixed_land"]]
    for taxon, src in prc_sources.items():
        if not src.is_file():
            raise FileNotFoundError(f"Missing final locked {taxon} PRC history: {src}")
        t = pd.read_csv(src, dtype={"NUTS_ID": str},
                        usecols=["NUTS_ID", "year", "prc_mean", "snow_status"])
        if not t.snow_status.eq("locked").all():
            raise ValueError(f"{taxon} PRC annual table lacks locked snow status")
        t = t[t.NUTS_ID.isin(selected)&t.year.isin(YEARS)].copy()
        if len(t) != 3*len(YEARS) or t.duplicated(["NUTS_ID", "year"]).any():
            raise ValueError(f"Selected {taxon} PRC annual trajectories incomplete")
        t["prc_mean"] = pd.to_numeric(t.prc_mean, errors="coerce")
        if not np.isfinite(t.prc_mean).all() or not t.prc_mean.between(0, 1).all():
            raise ValueError(f"Invalid final {taxon} annual PRC values")
        result = result.merge(
            t[["NUTS_ID", "year", "prc_mean"]].rename(
                columns={"prc_mean": f"{taxon}_prc"}
            ), on=["NUTS_ID", "year"], validate="one_to_one"
        )
    if len(result) != 3*len(YEARS):
        raise ValueError("Joining final annual sources changed selected regions")
    base = result[result.year.eq(2014)].set_index("NUTS_ID")
    for col in ("density_jan1_fixed_land", "spider_prc", "moth_prc"):
        result[f"change_{col}_from_2014"] = (
            result[col] - result.NUTS_ID.map(base[col])
        )
    # Validate 2014/2020 trajectories against original *frozen* endpoint source.
    for code, nid in ids.items():
        frozen = endpoints.loc[endpoints.NUTS_ID.eq(nid)].iloc[0]
        old = result.loc[result.NUTS_ID.eq(nid)&result.year.eq(2014)].iloc[0]
        final = result.loc[result.NUTS_ID.eq(nid)&result.year.eq(2020)].iloc[0]
        for year, rec in ((2014, old), (2020, final)):
            col = f"density_jan1_fixed_land_{year}"
            if not np.isclose(rec.density_jan1_fixed_land, frozen[col],
                              atol=1e-7, rtol=1e-9):
                raise ValueError(f"{nid} {year}: annual density != frozen endpoint")
        for taxon, endpoint_col in PRC.items():
            delta = final[f"change_{taxon}_prc_from_2014"]
            if not np.isclose(delta, frozen[endpoint_col], atol=2e-6, rtol=1e-5):
                raise ValueError(
                    f"{nid}: final annual {taxon} PRC Δ differs from frozen endpoint"
                )
    result["country"] = result.NUTS_ID.map({v:k for k,v in ids.items()})
    return result.sort_values(["country", "year"]).reset_index(drop=True)


def plot_4d(annual: pd.DataFrame, ids: dict[str, str], out: Path):
    order = ("FR", "TR", "NL")
    fig, axes = plt.subplots(3, 2, figsize=(figure_size("slide")[0], 7.85),
                             sharex=True, squeeze=False)
    for row, code in enumerate(order):
        g = annual[annual.NUTS_ID.eq(ids[code])].sort_values("year")
        left, right = axes[row]
        left.plot(g.year, g.change_density_jan1_fixed_land_from_2014,
                  color=SERIES_COLOURS["population"], linewidth=1.6,
                  marker="o", markersize=3.2)
        left.axhline(0, color="#666666", lw=.65)
        right.axhline(0, color="#666666", lw=.65)
        for taxon in ("spider", "moth"):
            right.plot(
                g.year, g[f"change_{taxon}_prc_from_2014"],
                color=SERIES_COLOURS[taxon], linewidth=1.6,
                marker="o", markersize=3.2
            )
        name = "Delft en Westland" if code == "NL" and ids[code] == "NL333" else ids[code]
        left.set_title(f"{COUNTRY_NAMES[code]} · {name}  |  Jan-1 density",
                       loc="left", fontsize=9)
        right.set_title(f"{COUNTRY_NAMES[code]} · {ids[code]}  |  PRC",
                        loc="left", fontsize=9)
        for ax in (left, right):
            ax.set_xticks(YEARS)
            ax.grid(axis="y", alpha=.17)
    axes[2,0].set_xlabel("Year")
    axes[2,1].set_xlabel("Year")
    fig.supylabel(LABEL_POP, x=.012, y=.55, fontsize=9)
     # Explicit per-axis labels avoid ambiguous dual-axis plotting.
    for ax in axes[:,0]:
        ax.set_ylabel("Δpeople km⁻²", fontsize=8)
    for ax in axes[:,1]:
        ax.set_ylabel("ΔPRC", fontsize=8)
    fig.legend(
        handles=[
            Line2D([0],[0],color=SERIES_COLOURS["population"],marker="o",label="Population density"),
            Line2D([0],[0],color=SERIES_COLOURS["spider"],marker="o",label="Spider PRC"),
            Line2D([0],[0],color=SERIES_COLOURS["moth"],marker="o",label="Moth PRC"),
        ],
        ncol=3, frameon=False, loc="upper right", fontsize=8,
        bbox_to_anchor=(.99, 1.015)
    )
    fig.suptitle("4D  |  Selected regional trajectories (2014 baseline = 0)",
                 x=.02, y=1.012, ha="left", fontsize=10)
    fig.subplots_adjust(left=.12, right=.98, top=.91, bottom=.075,
                        wspace=.38, hspace=.42)
    save_figure(fig, out)
    plt.close(fig)


def self_test():
    p = np.array([2., 3., -2., -3., 0., 1.])
    y = np.array([-1., 1., -1., 1., 2., 0.])
    assert categories(p, y, 0., 0.).tolist() == [
        "growth_prc_loss", "growth_prc_gain", "decline_prc_loss",
        "decline_prc_gain", "neutral", "neutral"
    ]
    strict = categories(np.array([4., 6.]), np.array([.0009, -.002]), 5., .001)
    assert strict.tolist() == ["neutral", "growth_prc_loss"]
    toy = pd.DataFrame({
        "NUTS_ID": ["FR001", "FR002", "TR001", "TR002", "NL333"],
        "country": ["FR", "FR", "TR", "TR", "NL"],
        POP: [10., 20., 15., 25., 8.],
        PRC["spider"]: [.01, .02, -.02, -.01, .01],
        PRC["moth"]: [.02, .04, -.04, -.02, .02],
    })
    labeled = classified(toy, 0., 0.)
    assert len(count_table(labeled)) == 10
    fr_id, fr_candidates = central_example(labeled, "FR", +1,
                                            pop_deadband=0., prc_deadband=0.)
    tr_id, tr_candidates = central_example(labeled, "TR", -1,
                                            pop_deadband=0., prc_deadband=0.)
    assert fr_id in {"FR001", "FR002"} and tr_id in {"TR001", "TR002"}
    assert fr_candidates == tr_candidates == 2
    assert all(count_table(labeled).groupby("taxon").n_nuts3.sum().eq(5))
    print("Figure 4 full sign partition, neutral rule and examples: PASS", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=ENDPOINTS)
    p.add_argument("--annual-pop", type=Path, default=ANNUAL_POP)
    p.add_argument("--spider-annual", type=Path, default=ANNUAL_PRC["spider"])
    p.add_argument("--moth-annual", type=Path, default=ANNUAL_PRC["moth"])
    p.add_argument("--out", type=Path, default=OUTPUT)
    p.add_argument("--pop-deadband", type=float, default=0.,
                   help="Neutral region when |Δdensity| <= this many people/km²")
    p.add_argument("--prc-deadband", type=float, default=0.,
                   help="Neutral region when |ΔPRC| <= this magnitude")
    p.add_argument("--fr-id", default=None,
                   help="Optional fixed French NUTS2016 NUTS3 for panel 4D")
    p.add_argument("--tr-id", default=None,
                   help="Optional fixed Turkish NUTS2016 NUTS3 for panel 4D")
    p.add_argument("--nl-id", default="NL333",
                   help="Optional Dutch case; default NL333")
    p.add_argument("--display-simplify-m", type=float, default=150.,
                   help="Vector map display-only EPSG:3035 simplify tolerance")
    p.add_argument("--expected-n", type=int, default=MATCHED_N,
                   help="Expected region count: 1273 (default) or verified 1452 UK extension")
    p.add_argument("--country-codes", nargs="+", default=list(COUNTRIES),
                   help="Countries to show in panel 4C (include UK for expanded cohort)")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    self_test()
    if args.self_test:
        return
    if not np.isfinite([args.pop_deadband,args.prc_deadband]).all() or (
        args.pop_deadband < 0 or args.prc_deadband < 0
    ):
        p.error("Classification deadbands must be finite and nonnegative")

    source = resolve(args.input)
    d = classified(frozen_endpoints(source, expected_n=args.expected_n),
                   args.pop_deadband, args.prc_deadband)
    counts = count_table(d)
    ids, chosen = example_ids(d, args)
    annual_sources = {
        "spider": resolve(args.spider_annual), "moth": resolve(args.moth_annual)
    }
    annual = trajectories(ids, d, resolve(args.annual_pop), annual_sources)
    geom, bounds = load_display_geometry(load_config(), args.display_simplify_m)
    out = resolve(args.out)
    out.mkdir(parents=True, exist_ok=True)
    apply_figure_style()
    # Source-aligned classifications are small audited derived tables.
    d[["NUTS_ID", "country", POP, *PRC.values(),
       "category_spider", "category_moth"]].to_csv(
        out/"fig4_nuts3_signed_category_audit.csv", index=False
    )
    counts.to_csv(out/"fig4A_category_counts.csv", index=False)
    annual.to_csv(out/"fig4D_selected_regional_annual_values.csv", index=False)
    plot_4a(counts, out/"fig4A_spider_moth_four_quadrant_counts", len(d))
    mapped = plot_4b(geom, bounds, d, out/"fig4B_europe_nuts3_sign_categories")
    if len(set(args.country_codes)) != len(args.country_codes):
        p.error("Duplicated --country-codes")
    by_country = plot_4c(d, tuple(args.country_codes),
                         out/"fig4C_country_category_proportions")
    by_country.to_csv(out/"fig4C_country_category_proportions.csv", index=False)
    plot_4d(annual, ids, out/"fig4D_fr_tr_nl333_regional_trajectories")
    manifest = {
        "source_frozen_matched_endpoints": str(source),
        "source_jan1_annual_density": str(resolve(args.annual_pop)),
        "source_final_taxon_annual_prc": {k:str(v) for k,v in annual_sources.items()},
        "cohort": "NUTS2016 Jan1 fixed land area",
        "years": [2014,2020],
        "matched_nuts3": len(d),
        "expected_nuts3_explicit": args.expected_n,
        "panel_4c_country_codes": list(args.country_codes),
        "source_prc_transformation": "locked native-pixel VIIRS->PRC before NUTS3 averaging",
        "pop_deadband_people_per_km2": args.pop_deadband,
        "prc_deadband": args.prc_deadband,
        "classification": "sign of ΔJan1 density vs sign of Δtaxon PRC",
        "neutral_definition": "either signed magnitude <= chosen neutral deadband",
        "n_neutral_by_taxon": {
            t:int(d[f"category_{t}"].eq("neutral").sum()) for t in PRC
        },
        "taxon_category_colours": COLORS,
        "map_display_extent_4326": list(EXTENT),
        "map_visual_simplification_epsg3035_metres": args.display_simplify_m,
        "mainland_polygon_counts_by_taxon": mapped,
        "selected_region_ids_for_D": chosen,
        "selection_note": "FR positive density and spider PRC gain; TR positive density and spider PRC loss, each closest to conditional bivariate medians; NL333 previously studied. Override IDs before final layout if desired.",
        "correlations_not_reestimated": True,
        "no_new_statistical_models": True,
        "source_data_unmodified": True,
        "outputs": "editable SVG and vector PDF",
    }
    (out/"fig4_provenance.json").write_text(
        json.dumps(manifest, indent=2)+"\n", encoding="utf-8"
    )
    print("\nFIGURE 4 OUTPUTS", flush=True)
    for artifact in sorted(out.glob("fig4*")):
        print(f"  {artifact} ({artifact.stat().st_size/1048576:.2f} MiB)", flush=True)
    print("Frozen matched regions:", len(d), flush=True)
    print("Figure 4D selected regions:", ids, flush=True)
    print("Classification neutral cases:",
          {t:int(d[f"category_{t}"].eq("neutral").sum()) for t in PRC},
          flush=True)


if __name__ == "__main__":
    main()
