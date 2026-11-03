#!/usr/bin/env python3
"""Publication Figure 2A–D: frozen NUTS2016 population/radiance comparison.

Read-only figure rendering from the EXISTING finished 1,273-region Jan1
2014–2020 country-analysis table. Do not re-extract VIIRS, reconstruct
population, redefine a cohort, or change any source model.

A: 2014 NUTS3 log10(Jan1 density) vs log10(corrected DNB radiance).
B: Matched 2014 and 2020 log-log clouds and descriptive linear trend lines.
C: All NUTS3 2014–2020 signed absolute density change vs signed Δ radiance.
D: Same Δpopulation vs Δradiance in existing selected-country facets:
   TR, DE, FR, IE, IT, NL.

SVG and vector PDF output; original scientific tables remain untouched.
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

from figure_style import apply_figure_style, figure_size, save_figure
from project_config import PROJECT_ROOT

ROOT = PROJECT_ROOT
INPUT = (ROOT / "data/processed/jan1_population_radiance_prc_country_nl333"
         / "matched_region_change_2014_2020.csv")
STATS = (ROOT / "data/processed/jan1_population_radiance_prc_country_nl333"
         / "within_country_correlations_2014_2020.csv")
OUTPUT = ROOT / "figures/publication/figure2_population_radiance"

COHORT_N = 1273
COUNTRIES = ("TR", "DE", "FR", "IE", "IT", "NL")
COUNTRY_LABELS = {
    "TR": "Türkiye", "DE": "Germany", "FR": "France",
    "IE": "Ireland", "IT": "Italy", "NL": "Netherlands",
    "UK": "United Kingdom",
}
DENSITY_COLS = {
    2014: "density_jan1_fixed_land_2014",
    2020: "density_jan1_fixed_land_2020",
}
RADIANCE_COLS = {2014: "radiance_2014", 2020: "radiance_2020"}
DX = "change_density_jan1_fixed_land_2014_2020"
DY = "radiance_delta_2014_2020"
YEAR_COLOURS = {2014: "#0072B2", 2020: "#D55E00"}
SCATTER_COLOUR = "#377BA6"
COUNTRY_COLOURS = {
    "TR": "#0072B2", "DE": "#E69F00", "FR": "#009E73",
    "IE": "#CC79A7", "IT": "#D55E00", "NL": "#7B61A8",
}
RAD_UNIT = r"nW cm$^{-2}$ sr$^{-1}$"
X_CHANGE = "Change in Jan-1 population density, 2014–2020 (people km⁻²)"
Y_CHANGE = "Change in calibrated VIIRS radiance, 2014–2020 (nW cm⁻² sr⁻¹)"
X_LOG = "log10(Jan-1 population density, people km⁻²)"
Y_LOG = "log10(calibrated VIIRS DNB radiance, nW cm⁻² sr⁻¹)"


def resolve(p: str | Path) -> Path:
    x = Path(p)
    return x if x.is_absolute() else ROOT / x


def load_frozen(path: Path, *, expected_n: int = COHORT_N) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(
            f"Frozen Jan-1 matched NUTS3 data missing: {path}. "
            "Do not regenerate or substitute an older Eurostat density cohort."
        )
    d = pd.read_csv(path, dtype={"NUTS_ID": str, "country": str})
    required = {"NUTS_ID", "country", DX, DY,
                *DENSITY_COLS.values(), *RADIANCE_COLS.values()}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"Frozen joined table lacks columns {sorted(missing)}")
    if len(d) != expected_n or not d.NUTS_ID.is_unique or d.NUTS_ID.isna().any():
        raise ValueError(f"Expected {expected_n} unique NUTS2016 regions; got {len(d)}")
    cols = sorted(required - {"NUTS_ID", "country"})
    for col in cols:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    if not np.isfinite(d[cols].to_numpy(float)).all():
        raise ValueError("Frozen NUTS3 source has missing or non-finite numeric values")
    if not d.country.str.fullmatch(r"[A-Z]{2}").all():
        raise ValueError("Unexpected country labels in frozen NUTS3 table")
    for y in DENSITY_COLS:
        if (d[DENSITY_COLS[y]] <= 0).any():
            raise ValueError(f"Non-positive reconstructed population density in {y}")
        if (d[RADIANCE_COLS[y]] < 0).any():
            raise ValueError(f"Negative calibrated radiance in {y}")
    if not np.allclose(
        d[DX], d[DENSITY_COLS[2020]] - d[DENSITY_COLS[2014]],
        atol=1e-9, rtol=1e-10
    ):
        raise ValueError("Frozen Δdensity is inconsistent with Jan1 endpoints")
    if not np.allclose(
        d[DY], d[RADIANCE_COLS[2020]] - d[RADIANCE_COLS[2014]],
        atol=1e-8, rtol=1e-8
    ):
        raise ValueError("Frozen Δradiance is inconsistent with VIIRS endpoints")
    return d.sort_values(["country", "NUTS_ID"]).reset_index(drop=True)


def matched_log_data(d: pd.DataFrame) -> pd.DataFrame:
    """Same finite, strictly positive NUTS3 set for BOTH log-log years."""
    positive = (d[RADIANCE_COLS[2014]] > 0) & (d[RADIANCE_COLS[2020]] > 0)
    eligible = d.loc[positive].copy()
    if len(eligible) < 1000:
        raise ValueError("Unexpectedly small positive-radiance log-log cohort")
    for year in (2014, 2020):
        eligible[f"log_density_{year}"] = np.log10(eligible[DENSITY_COLS[year]])
        eligible[f"log_radiance_{year}"] = np.log10(eligible[RADIANCE_COLS[year]])
    columns = [f"log_{metric}_{year}" for metric in ("density", "radiance")
               for year in (2014, 2020)]
    if not np.isfinite(eligible[columns].to_numpy(float)).all():
        raise ValueError("Non-finite matched logarithms")
    return eligible


def descriptive_fit(x, y, *, year: int) -> dict:
    """Visual OLS trend line (same descriptive approach as earlier plot scripts).

    The observational dataset remains frozen; this is a display line, not a
    new explanatory model or a replacement for prior country statistics.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) <= 0:
        raise ValueError(f"Insufficient variation for {year} visual trend line")
    slope, intercept = np.polyfit(x, y, 1)
    fitted = intercept + slope * x
    denom = np.sum((y - y.mean()) ** 2)
    r2 = (1 - np.sum((y - fitted) ** 2) / denom) if denom > 0 else np.nan
    return {
        "year": year, "n": len(x), "slope": float(slope),
        "intercept": float(intercept), "descriptive_r2": float(r2)
    }


def fit_y(fit: dict, x: np.ndarray) -> np.ndarray:
    return fit["intercept"] + fit["slope"] * x


def common_log_limits(log: pd.DataFrame) -> tuple[tuple[float, float], tuple[float, float]]:
    xx = np.concatenate([log["log_density_2014"], log["log_density_2020"]])
    yy = np.concatenate([log["log_radiance_2014"], log["log_radiance_2020"]])
    def lim(v):
        lo, hi = float(np.min(v)), float(np.max(v))
        pad = .045 * max(hi - lo, .1)
        return lo-pad, hi+pad
    return lim(xx), lim(yy)


def plot_2a(log: pd.DataFrame, fit2014: dict, limits, out: Path):
    fig, ax = plt.subplots(figsize=figure_size("double", .70))
    x = log["log_density_2014"].to_numpy(float)
    y = log["log_radiance_2014"].to_numpy(float)
    ax.scatter(x, y, s=9, alpha=.45, color=YEAR_COLOURS[2014],
               linewidths=0, label=f"2014 NUTS3 (n={len(log):,})")
    line_x = np.linspace(limits[0][0], limits[0][1], 200)
    ax.plot(line_x, fit_y(fit2014, line_x), lw=1.8, color="#222222",
            label=f"OLS display line (slope={fit2014['slope']:.2f})")
    ax.set(xlabel=X_LOG, ylabel=Y_LOG, title="2A  |  Spatial relationship, 2014")
    ax.set_xlim(limits[0]); ax.set_ylim(limits[1])
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=.14)
    fig.tight_layout()
    save_figure(fig, out)
    plt.close(fig)


def plot_2b(log: pd.DataFrame, fits: dict, limits, out: Path):
    fig, ax = plt.subplots(figsize=figure_size("double", .70))
    line_x = np.linspace(limits[0][0], limits[0][1], 200)
    for year in (2014, 2020):
        x = log[f"log_density_{year}"].to_numpy(float)
        y = log[f"log_radiance_{year}"].to_numpy(float)
        ax.scatter(x, y, s=8, alpha=.24, color=YEAR_COLOURS[year],
                   linewidths=0)
        fit = fits[year]
        ax.plot(line_x, fit_y(fit, line_x), lw=1.8,
                color=YEAR_COLOURS[year],
                label=f"{year} (slope={fit['slope']:.2f}; n={len(log):,})")
    ax.set(xlabel=X_LOG, ylabel=Y_LOG,
           title="2B  |  Spatial relationships, 2014 and 2020")
    ax.set_xlim(limits[0]); ax.set_ylim(limits[1])
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=.14)
    fig.tight_layout()
    save_figure(fig, out)
    plt.close(fig)


def frozen_rho(csv: Path, d: pd.DataFrame):
    """OPTIONAL only: display already-frozen Spearman statistics.

    Do not calculate a new population/radiance correlation in this figure
    script. If the original frozen CSV is not available, omit that label.
    """
    if not csv.is_file():
        print(f"NOTE: Frozen within-country table missing; omit ρ annotations: {csv}")
        return {}
    a = pd.read_csv(csv, dtype={"country": str, "outcome": str})
    needed = {"country", "outcome", "n", "spearman_rho"}
    if not needed.issubset(a):
        raise ValueError(f"Frozen correlation table missing {needed-set(a)}")
    r = a[a.outcome.eq("radiance")]
    if r.country.duplicated().any():
        raise ValueError("Repeated country radiance outcomes in frozen statistics")
    result = {}
    for _, row in r.iterrows():
        country = str(row.country)
        n = len(d) if country == "ALL" else int(d.country.eq(country).sum())
        if n != int(row.n):
            raise ValueError(
                f"Frozen country stats cohort mismatch: {country} {row.n} vs {n}"
            )
        result[country] = float(row.spearman_rho)
    return result


def signed_limits(x: np.ndarray, padding: float = .1) -> tuple[float, float]:
    vmax = float(np.max(np.abs(np.asarray(x, float))))
    vmax = vmax * (1 + padding)
    return -max(vmax, 1e-8), max(vmax, 1e-8)


def setup_change_axes(ax, *, xlin: float, ylin: float, limits):
    ax.set_xscale("symlog", linthresh=xlin)
    ax.set_yscale("symlog", linthresh=ylin)
    ax.set_xlim(limits[0])
    ax.set_ylim(limits[1])
    ax.axhline(0, color="#777777", lw=.65, zorder=0)
    ax.axvline(0, color="#777777", lw=.65, zorder=0)
    ax.grid(alpha=.13, which="major")


def plot_2c(d: pd.DataFrame, rho: dict, *, xlin: float, ylin: float, limits, out: Path):
    fig, ax = plt.subplots(figsize=figure_size("double", .70))
    ax.scatter(d[DX], d[DY], s=10, linewidths=0, alpha=.42,
               color=SCATTER_COLOUR)
    setup_change_axes(ax, xlin=xlin, ylin=ylin, limits=limits)
    ax.set(xlabel=X_CHANGE, ylabel=Y_CHANGE,
           title="2C  |  Regional changes, 2014–2020")
    subtitle = f"n={len(d):,} NUTS3; signed-log axes"
    if "ALL" in rho and np.isfinite(rho["ALL"]):
        subtitle += f"; frozen Spearman ρ={rho['ALL']:+.2f}"
    ax.text(.015, .98, subtitle, transform=ax.transAxes, va="top",
            fontsize=8, color="#444444")
    fig.tight_layout()
    save_figure(fig, out)
    plt.close(fig)


def plot_2d(d: pd.DataFrame, rho: dict, countries: tuple[str, ...],
            *, xlin: float, ylin: float, limits, out: Path):
    ncol = 2
    nrow = int(np.ceil(len(countries)/ncol))
    # Extra height ensures readable 3 × 2 facets on a single editable sheet.
    fig, axs = plt.subplots(nrow, ncol,
                           figsize=(figure_size("slide")[0], 3.18*nrow),
                           sharex=True, sharey=True, squeeze=False)
    for i, code in enumerate(countries):
        ax = axs.flat[i]
        sub = d[d.country.eq(code)]
        if sub.empty:
            raise ValueError(f"Country {code} is absent from frozen data")
        ax.scatter(sub[DX], sub[DY], s=12, alpha=.62,
                   linewidths=0, color=COUNTRY_COLOURS.get(code, SCATTER_COLOUR))
        # NL333 is the established Dutch outlier; identify but do not remove it.
        if code == "NL" and "NL333" in set(sub.NUTS_ID):
            focal = sub.loc[sub.NUTS_ID.eq("NL333")].iloc[0]
            ax.scatter([focal[DX]], [focal[DY]], marker="*",
                       s=105, c="#222222", zorder=6)
            ax.annotate("NL333", (focal[DX], focal[DY]),
                        xytext=(4, 4), textcoords="offset points", fontsize=7)
        setup_change_axes(ax, xlin=xlin, ylin=ylin, limits=limits)
        title = f"{COUNTRY_LABELS.get(code, code)} ({code}); n={len(sub)}"
        if code in rho and np.isfinite(rho[code]):
            title += f"; ρ={rho[code]:+.2f}"
        ax.set_title(title, fontsize=9, pad=8, loc="left")
    for ax in list(axs.flat)[len(countries):]:
        ax.set_visible(False)
    fig.suptitle("2D  |  Country-specific NUTS3 changes, 2014–2020", y=.997, fontsize=11)
    fig.supxlabel(X_CHANGE, fontsize=9, y=.025)
    fig.supylabel(Y_CHANGE, fontsize=9, x=.01)
    fig.subplots_adjust(top=.935, bottom=.08, left=.12, right=.98,
                        hspace=.28, wspace=.22)
    save_figure(fig, out)
    plt.close(fig)


def self_test():
    data = pd.DataFrame({
        "NUTS_ID": ["X", "Y", "Z"], "country": ["FR", "DE", "NL"],
        DENSITY_COLS[2014]: [10., 100., 300.],
        DENSITY_COLS[2020]: [12., 101., 310.],
        RADIANCE_COLS[2014]: [.1, 3., 8.],
        RADIANCE_COLS[2020]: [.2, 2., 9.],
        DX: [2., 1., 10.], DY: [.1, -1., 1.],
    })
    validated = load_frozen_sample(data)
    assert len(validated) == 3
    fit = descriptive_fit(np.array([1., 2., 3.]),
                          np.array([2., 4., 6.]), year=2014)
    assert np.isclose(fit["slope"], 2) and np.isclose(fit["intercept"], 0)
    assert signed_limits(np.array([-2., 3.])) == (-3.3, 3.3) or (
        np.allclose(signed_limits(np.array([-2., 3.])), (-3.3, 3.3))
    )
    print("Figure 2 paired-endpoint, units and visual-line self-test: PASS", flush=True)


def load_frozen_sample(data: pd.DataFrame) -> pd.DataFrame:
    """Tiny independent unit test without accessing user-specific data."""
    x = data.copy()
    assert np.allclose(x[DX], x[DENSITY_COLS[2020]]-x[DENSITY_COLS[2014]])
    assert np.allclose(x[DY], x[RADIANCE_COLS[2020]]-x[RADIANCE_COLS[2014]])
    assert x.NUTS_ID.is_unique
    return x


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=INPUT)
    p.add_argument("--country-stats", type=Path, default=STATS)
    p.add_argument("--out", type=Path, default=OUTPUT)
    p.add_argument("--country-codes", nargs="+", default=list(COUNTRIES))
    p.add_argument("--density-linthresh", type=float, default=5.,
                   help="Signed-log change x-axis linear region, people/km²")
    p.add_argument("--radiance-linthresh", type=float, default=.1,
                   help="Signed-log change y-axis linear region, nW cm^-2 sr^-1")
    p.add_argument("--expected-n", type=int, default=COHORT_N,
                   help="Expected unique region count (1273 frozen; 1452 explicit UK extension)")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    self_test()
    if args.self_test:
        return
    if args.density_linthresh <= 0 or args.radiance_linthresh <= 0:
        p.error("Signed-log thresholds must be positive")
    countries = tuple(args.country_codes)
    if not countries or len(countries) != len(set(countries)):
        p.error("Provide distinct country codes")
    source, stats_file, out = resolve(args.input), resolve(args.country_stats), resolve(args.out)
    d = load_frozen(source, expected_n=args.expected_n)
    log = matched_log_data(d)
    rho = frozen_rho(stats_file, d)
    fits = {
        y: descriptive_fit(
            log[f"log_density_{y}"], log[f"log_radiance_{y}"], year=y
        ) for y in (2014, 2020)
    }
    log_limits = common_log_limits(log)
    change_limits = (signed_limits(d[DX].to_numpy(float)),
                     signed_limits(d[DY].to_numpy(float)))
    apply_figure_style()
    out.mkdir(parents=True, exist_ok=True)
    plot_2a(log, fits[2014], log_limits, out/"fig2A_spatial_population_radiance_2014")
    plot_2b(log, fits, log_limits, out/"fig2B_spatial_relation_2014_2020")
    plot_2c(d, rho, xlin=args.density_linthresh, ylin=args.radiance_linthresh,
            limits=change_limits, out=out/"fig2C_population_vs_radiance_change_2014_2020")
    plot_2d(d, rho, countries, xlin=args.density_linthresh,
            ylin=args.radiance_linthresh, limits=change_limits,
            out=out/"fig2D_country_population_radiance_change_2014_2020")
    # Only plotting regressions are newly evaluated; all country ρ labels are
    # taken from the pre-existing frozen output CSV.
    pd.DataFrame(list(fits.values())).to_csv(out/"fig2_descriptive_plot_line_parameters.csv",
                                               index=False)
    manifest = {
        "frozen_input": str(source),
        "frozen_correlations": str(stats_file) if stats_file.is_file() else None,
        "n_nuts3_all_signed_changes": len(d),
        "expected_nuts3_explicit": args.expected_n,
        "n_nuts3_matched_log_positive_radiance": len(log),
        "n_nuts3_not_log_plotted": len(d)-len(log),
        "years": [2014, 2020],
        "population_definition": "Eurostat 1 Jan population / fixed NUTS2016 NUTS3 land area",
        "radiance_definition": "calibrated annual VIIRS DNB NUTS2016 regional mean",
        "visual_line_method": "ordinary least squares on paired positive log-log set; for graphic description only",
        "figure_C_D_axes": "symmetric logarithmic, unchanged physical signed Δ variables",
        "density_linthresh_people_km2": args.density_linthresh,
        "radiance_linthresh_nw_cm2_sr": args.radiance_linthresh,
        "country_codes_D": list(countries),
        "D_country_spearman_source": "frozen within-country table; not recalculated",
        "country_spearman_radiance": {k: v for k, v in rho.items()
                                     if k == "ALL" or k in countries},
        "unmodified_scientific_inputs": True,
        "vector_formats": ["svg", "pdf"],
    }
    (out/"fig2_provenance.json").write_text(json.dumps(manifest, indent=2)+"\n",
                                            encoding="utf-8")
    print("\nFIGURE 2 OUTPUTS", flush=True)
    for name in sorted(out.glob("fig2*")):
        print(f"  {name}  ({name.stat().st_size/1048576:.2f} MiB)", flush=True)
    print(f"Frozen cohort: n={len(d):,}; matched positive log-log: n={len(log):,}", flush=True)


if __name__ == "__main__":
    main()
