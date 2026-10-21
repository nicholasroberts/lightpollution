#!/usr/bin/env python3
"""Matched 2014–2020 Jan1-density, VIIRS and spider/moth PRC country comparisons.

Reuses the SAME 1,273 fixed-NUTS2016 regions from the existing four-panel map.
Jan1 population/fixed land area, locked biological response and VIIRS history are
READ ONLY. No PRC refits or new population harmonisation. The primary analyses
are continuous signed changes, not categories chosen to achieve a hypothesis.
NL333 is compared with the rest of NL and actual neighbouring NUTS3 polygons.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from project_config import PROJECT_ROOT as ROOT, load_config
from historical_prc_pixel_nuts3 import load_nuts
from plot_jan1_population_prc_change_symlog import (
    INPUT, RADIANCE_SOURCE, read_matched_endpoints, add_radiance_endpoints,
    robust_symmetric_limit
)

YEAR = 2020
OUT = ROOT / "data/processed/jan1_population_radiance_prc_country_nl333"
FIG = ROOT / "figures/jan1_population_radiance_prc_country_nl333"
METRICS = {
    "radiance": ("radiance_delta_2014_2020", "VIIRS radiance change (nW cm⁻² sr⁻¹)"),
    "spider": ("cumulative_delta_spider_prc_2020", "Spider PRC change"),
    "moth": ("cumulative_delta_moth_prc_2020", "Moth PRC change"),
}
POP = "change_density_jan1_fixed_land_2014_2020"
DISPLAY = ("TR", "DE", "FR", "IE", "IT", "NL")


def matched_data(endpoints: Path, radiance: Path) -> pd.DataFrame:
    """Share the four-panel script's validated population and radiance join."""
    x = add_radiance_endpoints(read_matched_endpoints(endpoints), radiance)
    need = [POP, "fixed_land_area_km2", "radiance_2014",
            *[item[0] for item in METRICS.values()]]
    if not np.isfinite(x[need].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite matched analysis values")
    if (x.fixed_land_area_km2 <= 0).any():
        raise ValueError("Non-positive fixed regional land area")
    if x.NUTS_ID.nunique() != 1273 or not x.NUTS_ID.is_unique:
        raise ValueError("Changed fixed reconstructed cohort")
    if not x.country.str.len().eq(2).all():
        raise ValueError("Expected ISO-like two-letter country identifiers")
    if "NL333" not in set(x.NUTS_ID):
        raise ValueError("NL333 is absent from matched Jan1 cohort")
    return x


def stats(x: pd.Series, y: pd.Series) -> dict:
    """Descriptive associations, without misleading independent-region p-values."""
    a, b = np.asarray(x, float), np.asarray(y, float)
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if len(a) < 4 or np.unique(a).size < 2 or np.unique(b).size < 2:
        return {"n": len(a), "spearman_rho": np.nan, "pearson_r": np.nan}
    return {"n": len(a), "spearman_rho": float(spearmanr(a, b).statistic),
            "pearson_r": float(pearsonr(a, b).statistic)}


def within_countries(x: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, sub in [("ALL", x), *list(x.groupby("country"))]:
        for label, (col, _) in METRICS.items():
            rows.append({"country": name, "outcome": label,
                         **stats(sub[POP], sub[col]),
                         "n_growth_over_5": int((sub[POP] > 5).sum()),
                         "n_decline_below_minus_5": int((sub[POP] < -5).sum()),
                         "median_density_change": float(sub[POP].median()),
                         "median_outcome_change": float(sub[col].median())})
    return pd.DataFrame(rows).sort_values(["outcome", "country"])


def country_aggregates(x: pd.DataFrame) -> pd.DataFrame:
    """Equal-country analysis uses one *area-weighted* endpoint delta per country.

    Fixed land area weights are applied consistently to all NUTS3 values.
    Country population-density delta is equivalently total Jan1 population
    change divided by total fixed covered land area. PRC/radiance are averages
    of region-level quantities (native VIIRS pixels were transformed BEFORE
    PRC aggregation), not PRC evaluated at mean country radiance.
    """
    rows = []
    for country, sub in x.groupby("country"):
        w = sub.fixed_land_area_km2.to_numpy(dtype=float)
        row = {"country": country, "n_regions": len(sub),
               "matched_fixed_land_km2": float(w.sum()),
               "population_density_change_people_per_km2":
                   float(np.average(sub[POP], weights=w)),
               "region_median_population_change": float(sub[POP].median())}
        for label, (col, _) in METRICS.items():
            row[f"{label}_change_area_weighted"] = float(np.average(sub[col], weights=w))
            row[f"{label}_change_region_median"] = float(sub[col].median())
        rows.append(row)
    return pd.DataFrame(rows).sort_values("country")


def between_countries(countries: pd.DataFrame) -> pd.DataFrame:
    """Each country is one observation, not each of its NUTS3 subdivisions."""
    rows = []
    predictor = "population_density_change_people_per_km2"
    for min_regions in (1, 5, 10):
        eligible = countries[countries.n_regions >= min_regions]
        for label in METRICS:
            col = f"{label}_change_area_weighted"
            rows.append({"n_min_regions": min_regions, "outcome": label,
                         "aggregation": "fixed-land-area-weighted national changes; unweighted country-level correlation",
                         **stats(eligible[predictor], eligible[col])})
    return pd.DataFrame(rows)


def nl333_influence(x: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    nl = x[x.country.eq("NL")].copy()
    if len(nl) < 5:
        raise ValueError("Netherlands has insufficient matched NUTS3 regions")
    focal = nl[nl.NUTS_ID.eq("NL333")].iloc[0]
    others = nl[nl.NUTS_ID.ne("NL333")]
    rows = []
    for label, (col, _) in METRICS.items():
        base = stats(nl[POP], nl[col])
        reduced = stats(others[POP], others[col])
        absolute_rank = int(nl[col].abs().rank(method="min", ascending=False)
                            [nl.NUTS_ID.eq("NL333")].iloc[0])
        rows.append({"outcome": label, "n_nl": len(nl),
                     "nl333_density_change": float(focal[POP]),
                     "nl333_outcome_change": float(focal[col]),
                     "nl333_rank_absolute_outcome_in_nl": absolute_rank,
                     "nl333_baseline_radiance": float(focal.radiance_2014),
                     "rho_with_nl333": base["spearman_rho"],
                     "rho_without_nl333": reduced["spearman_rho"],
                     "pearson_with_nl333": base["pearson_r"],
                     "pearson_without_nl333": reduced["pearson_r"]})
    return pd.DataFrame(rows), nl


def local_geography(nl: pd.DataFrame, nuts) -> tuple[pd.DataFrame, object]:
    """Distances from ACTUAL NUTS2016 boundaries in EPSG:3035 (metres).

    Returns every matched Dutch region, classified by touching NL333 or
    shortest polygon-boundary distance. This is NUTS3-neighbourhood analysis,
    not a pixel/municipality-level attribution to horticultural lamps.
    """
    geo = nuts[["NUTS_ID", "geometry"]].copy()
    if geo.NUTS_ID.duplicated().any():
        raise ValueError("NUTS geometry contains duplicated region IDs")
    g = geo.merge(nl, on="NUTS_ID", how="inner", validate="one_to_one").to_crs("EPSG:3035")
    if len(g) != len(nl) or not g.NUTS_ID.eq("NL333").any():
        raise ValueError("Dutch geographic merge incomplete")
    anchor = g.loc[g.NUTS_ID.eq("NL333"), "geometry"].iloc[0]
    centre = anchor.centroid
    g["distance_boundary_km"] = g.geometry.distance(anchor) / 1000.
    g["distance_centroid_km"] = g.geometry.centroid.distance(centre) / 1000.
    g["touches_nl333"] = g.geometry.touches(anchor)
    g["within_25km_boundary"] = g.distance_boundary_km.le(25.)
    g["within_50km_boundary"] = g.distance_boundary_km.le(50.)
    g["is_nl333"] = g.NUTS_ID.eq("NL333")
    g = g.sort_values(["distance_boundary_km", "distance_centroid_km"])
    return pd.DataFrame(g.drop(columns="geometry")), g


def nl_local_correlations(local: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, sub in [
        ("NL all", local),
        ("NL without NL333", local[~local.is_nl333]),
        ("NL333 + polygon touching", local[local.touches_nl333 | local.is_nl333]),
        ("NL333 + within 25km of boundary", local[local.within_25km_boundary]),
        ("NL333 + within 50km of boundary", local[local.within_50km_boundary]),
    ]:
        for label, (col, _) in METRICS.items():
            rows.append({"group": group, "outcome": label, **stats(sub[POP], sub[col])})
    return pd.DataFrame(rows)


def scatter_figure(x: pd.DataFrame, out: Path) -> None:
    fig, axs = plt.subplots(len(DISPLAY), 3, figsize=(15, 20), layout="constrained")
    for row, country in enumerate(DISPLAY):
        sub = x[x.country.eq(country)]
        for col, (label, (y, y_label)) in enumerate(METRICS.items()):
            ax = axs[row, col]
            ax.axhline(0, lw=.6, color="gray")
            ax.axvline(0, lw=.6, color="gray")
            ax.scatter(sub[POP], sub[y], s=17, alpha=.7, color="steelblue")
            if country == "NL":
                z = sub[sub.NUTS_ID.eq("NL333")]
                ax.scatter(z[POP], z[y], marker="*", s=140,
                           color="darkorange", edgecolor="black", zorder=5)
                ax.annotate("NL333", (float(z[POP].iloc[0]), float(z[y].iloc[0])),
                            xytext=(4, 5), textcoords="offset points", fontsize=8)
            rho = stats(sub[POP], sub[y])["spearman_rho"]
            ax.set_title(f"{country} | n={len(sub)} | Spearman rho={rho:.2f}")
            if row == len(DISPLAY) - 1:
                ax.set_xlabel("2014–2020 Δ Jan1 density (people/km²)")
            if col == 0:
                ax.set_ylabel("VIIRS Δ (nW cm⁻² sr⁻¹)")
            elif col == 1:
                ax.set_ylabel("Spider ΔPRC")
            else:
                ax.set_ylabel("Moth ΔPRC")
    fig.suptitle("Observed Jan1 density change vs corrected VIIRS and locked pixel PRC (2014–2020)", fontsize=14)
    fig.savefig(out / "selected_country_population_vs_radiance_spider_moth.pdf",
                bbox_inches="tight")
    plt.close(fig)


def between_figure(countries: pd.DataFrame, out: Path) -> None:
    fig, axs = plt.subplots(1, 3, figsize=(16, 5), layout="constrained")
    xcol = "population_density_change_people_per_km2"
    for ax, (label, (_, y_label)) in zip(axs, METRICS.items()):
        ycol = f"{label}_change_area_weighted"
        ax.axhline(0, color="gray", lw=.7)
        ax.axvline(0, color="gray", lw=.7)
        ax.scatter(countries[xcol], countries[ycol], s=25, color="steelblue", alpha=.7)
        for _, row in countries[countries.country.isin(DISPLAY)].iterrows():
            ax.annotate(row.country, (row[xcol], row[ycol]), xytext=(3, 4),
                        textcoords="offset points", fontsize=9)
        result = stats(countries[xcol], countries[ycol])
        ax.set_title(f"Between countries: {label}\nn={result['n']}, Spearman rho={result['spearman_rho']:.2f}")
        ax.set_xlabel("Country Δ Jan1 density (people/km²)")
        ax.set_ylabel(y_label)
    fig.suptitle("Equal weight per country; each axis value is fixed-land-weighted over the matched region cohort")
    fig.savefig(out / "between_country_population_vs_radiance_spider_moth.pdf",
                bbox_inches="tight")
    plt.close(fig)


def nl_focus_map(geo, out: Path) -> None:
    g = geo.copy()
    anchor = g[g.is_nl333].iloc[0].geometry
    cx, cy = anchor.centroid.x, anchor.centroid.y
    columns = [(POP, "Jan1 population Δ", "symlog", 1.),
               (METRICS["radiance"][0], "Corrected VIIRS Δ", "symlog", .1),
               (METRICS["spider"][0], "Spider ΔPRC", "linear", 1.),
               (METRICS["moth"][0], "Moth ΔPRC", "linear", 1.)]
    fig, axs = plt.subplots(1, 4, figsize=(19, 5.7), layout="constrained")
    for ax, (col, label, mode, threshold) in zip(axs, columns):
        limit = robust_symmetric_limit(g[col], 1e-10)
        norm = (SymLogNorm(linthresh=threshold, vmin=-limit, vmax=limit)
                if mode == "symlog" else
                TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
        g.plot(column=col, ax=ax, cmap="RdBu_r", norm=norm,
               edgecolor="#888888", linewidth=.25)
        g.loc[g.is_nl333].boundary.plot(ax=ax, color="black", linewidth=2)
        g.loc[g.touches_nl333].boundary.plot(ax=ax, color="orange", linewidth=1)
        ax.set(xlim=(cx-80000, cx+80000), ylim=(cy-80000, cy+80000),
               title=f"{label} | {mode}")
        ax.set_aspect("equal")
        ax.set_axis_off()
        sm = ScalarMappable(norm=norm, cmap="RdBu_r")
        sm.set_array([])
        fig.colorbar(sm, ax=ax, shrink=.65, extend="both")
    fig.suptitle("NL333 Delft en Westland (black) and boundary-sharing NUTS3 neighbours (orange), 2014–2020")
    fig.savefig(out / "nl333_four_panel_local_80km.pdf", bbox_inches="tight")
    fig.savefig(out / "nl333_four_panel_local_80km.png", bbox_inches="tight", dpi=200)
    plt.close(fig)


def synthetic_test() -> None:
    x = pd.Series([0., 1., 2., 3., 4.])
    y = pd.Series([0., 1., 2., 3., 4.])
    assert np.isclose(stats(x, y)["spearman_rho"], 1.)
    assert np.isnan(stats(x.iloc[:2], y.iloc[:2])["spearman_rho"])
    toy = pd.DataFrame({
        "NUTS_ID": ["NL333", "NL331", "NL332", "FR001", "FR002"],
        "country": ["NL", "NL", "NL", "FR", "FR"],
        "fixed_land_area_km2": [1., 1., 2., 2., 3.],
        POP: [1., 2., 3., -1., -2.],
        "radiance_2014": [10., 2., 1., 1., 2.],
        METRICS["radiance"][0]: [-2., 1., 2., -1., -2.],
        METRICS["spider"][0]: [.1, -.1, -.2, .1, .2],
        METRICS["moth"][0]: [.2, -.2, -.3, .1, .2],
    })
    c = country_aggregates(toy)
    nl = c[c.country.eq("NL")].iloc[0]
    assert np.isclose(nl.population_density_change_people_per_km2, 2.25)
    assert np.isclose(nl.radiance_change_area_weighted, .75)
    assert len(within_countries(toy)) == 9
    assert toy[toy.country.eq("NL")].NUTS_ID.nunique() == 3
    assert len(between_countries(c)) == 9
    assert np.isfinite(stats(toy[POP], toy[METRICS["radiance"][0]])["spearman_rho"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--radiance", type=Path, default=RADIANCE_SOURCE)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--figures", type=Path, default=FIG)
    args = parser.parse_args()
    synthetic_test()
    if args.self_test:
        print("Matched Jan1 country/NL333 synthetic checks: PASS")
        return
    x = matched_data(args.input, args.radiance)
    out, figures = args.out, args.figures
    out.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    within = within_countries(x)
    countries = country_aggregates(x)
    between = between_countries(countries)
    influence, nl = nl333_influence(x)
    nuts = load_nuts(load_config())
    local, local_geometry = local_geography(nl, nuts)
    local_corr = nl_local_correlations(local)
    x.to_csv(out / "matched_region_change_2014_2020.csv", index=False)
    within.to_csv(out / "within_country_correlations_2014_2020.csv", index=False)
    countries.to_csv(out / "country_area_weighted_changes_2014_2020.csv", index=False)
    between.to_csv(out / "between_country_correlations_2014_2020.csv", index=False)
    influence.to_csv(out / "nl333_influence_2014_2020.csv", index=False)
    local.to_csv(out / "nl333_neighbourhood_nuts3_2014_2020.csv", index=False)
    local_corr.to_csv(out / "nl333_neighbourhood_correlations_2014_2020.csv", index=False)
    scatter_figure(x, figures)
    between_figure(countries, figures)
    nl_focus_map(local_geometry, figures)
    manifest = {
        "years": [2014, 2020],
        "population": str(args.input), "radiance": str(args.radiance),
        "population_definition": "Jan1 population over region-specific FIXED land area",
        "matched_regions": len(x), "matched_countries": int(x.country.nunique()),
        "nl_regions": len(nl), "analysis_lock": "Input datasets treated as fixed; NOT modified",
        "comparison_unit_within": "NUTS2016 NUTS3 region",
        "comparison_unit_between": "country (equal country weight after land-area-weighted regional aggregation)",
        "nl333_location": "NUTS2016 Delft en Westland; distances among Dutch NUTS3 polygons in EPSG:3035",
        "geographic_neighbourhoods": "polygon touches / up to 25km or 50km shortest boundary distances",
        "warnings": [
            "Descriptive observational ranks cannot demonstrate causation or absence of effects",
            "Nominal independent-NUTS3 p-values intentionally omitted due to spatial dependence",
            "No 2014–2024 NL333 radiance claims transferred to the 2014–2020 interval",
            "Radiance and pixel-transformed PRC share VIIRS sensor data: not independent biological replication",
            "Country aggregate PRC uses consistent fixed land-area weighting of existing NUTS3 pixel-derived means",
            "Regional Jan1 demographic support may differ from VIIRS valid-pixel support",
            "NL333 geographical neighbours are not the same as individual municipal or greenhouse footprints"
        ]
    }
    (out / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Country comparisons (2014–2020, matched Jan1 cohort):")
    print(within[within.country.isin(["ALL", *DISPLAY])]
          [["country", "outcome", "n", "spearman_rho"]].to_string(index=False))
    print("\nNL333 influence:")
    print(influence.to_string(index=False))
    print(f"\nCSV: {out}\nFigures: {figures}")


if __name__ == "__main__":
    main()
