#!/usr/bin/env python3
"""Two publication-ready, intervention-annotated annual + monthly VIIRS figures.

A: frozen annual NUTS2016 VNL V2 mean 2014–2024 (verified 2021 snow patch).
B: native EOG VNL V1 monthly MUNICIPAL mean, 2020–2024 (Coventry 2022–24).
C: same-calendar-month MUNICIPAL VNL V1 anomaly for intervention year.

Annual NUTS3 and municipal monthly VNL products are displayed separately;
their different spatial footprints, products and compositing MUST NOT be
combined numerically. Event dates from independent city records, not inferred
by choosing visually convenient points.

   python src/plot_coventry_echirolles_annual_monthly_attribution.py --self-test
   python src/plot_coventry_echirolles_annual_monthly_attribution.py

If a validated annual bridge CSV does not exist, the plotter generates the
read-only annual table from the frozen parquet, endpoint and verified *existing*
2021 patch via export_coventry_echirolles_nuts3_bridge. No source changes.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from export_coventry_echirolles_nuts3_bridge import (
    ANNUAL, ENDPOINT, SNOW_CACHE, annual_series, endpoint_series,
    strict_snow_cache, REGIONS,
)
from local_viirs_five_case_analysis import ROOT

OUT = ROOT / "figures/local_lighting_interventions/nuts3_event_attribution"
DATA_OUT = ROOT / "data/processed/local_lighting_interventions/nuts3_event_attribution"
MONTHLY = (ROOT / "data/processed/local_lighting_interventions/five_case_monthly")
YEAR_FIRST, YEAR_LAST = 2014, 2024
RUNIT = "nW cm$^{-2}$ sr$^{-1}$"


@dataclass(frozen=True)
class Setup:
    key: str
    region: str
    title: str
    outcome_year: int
    references: tuple[int, ...]
    comparison_months: tuple[int, ...]
    first_month: str
    event_window: tuple[str, str]
    event_label: str
    dates: tuple[tuple[str, str, str], ...]
    led_window: tuple[str, str] | None
    notes: str


CASES = (
    Setup(
        key="coventry", region="UKG33", title="Coventry",
        outcome_year=2024, references=(2022, 2023),
        comparison_months=(8, 9, 10, 11, 12), first_month="2022-01",
        event_window=("2024-05-01", "2024-08-01"),
        event_label="Part-night switch-off rollout: May–July 2024",
        dates=(
            ("2024-02-20", "20 Feb: approval", "approval"),
            ("2024-05-01", "May: programming began", "rollout"),
            ("2024-08-01", "Aug: post-rollout", "post"),
        ),
        led_window=None,
        notes=("Approximately 60% of municipal streetlights eligible; "
               "earlier LEDs and dimming existed. Major LED programme "
               "started Nov 2025, OUTSIDE both series."),
    ),
    Setup(
        key="echirolles", region="FRK24", title="Échirolles / Isère",
        outcome_year=2022, references=(2020, 2021),
        comparison_months=(10, 11, 12), first_month="2020-01",
        event_window=("2022-07-01", "2023-01-01"),
        event_label="Phased switch-offs: Jul–Dec 2022",
        dates=(
            ("2022-07-29", "29 Jul: first sites", "rollout"),
            ("2022-08-01", "Aug: parks/buildings", "rollout"),
            ("2022-10-01", "Oct: road phase starts", "rollout"),
            ("2023-01-01", "Jan: road phase ended", "post"),
        ),
        led_window=("2023-06-01", "2025-01-01"),
        notes=("Échirolles is a subset of FRK24 Isère. LED upgrades "
               "increased from summer 2023; do not attribute Isère's "
               "entire decadal decline to this one city."),
    ),
)


def validated_annual():
    """Rebuild in MEMORY only; source Figure1 endpoint/snow lock is verified.

    No source files are written or overwritten. Reuses the existing exporter
    validation rather than trusting previously generated arbitrary CSVs.
    """
    if not ANNUAL.is_file() or not ENDPOINT.is_file():
        raise FileNotFoundError(
            f"Frozen source missing. Check {ANNUAL} and {ENDPOINT}; "
            "do not substitute older uncontrolled NUTS3 time series.")
    snow = strict_snow_cache()  # Fail closed on missing/stale source audit.
    raw = pd.read_parquet(
        ANNUAL, columns=["NUTS_ID", "year", "radiance_mean_corrected"])
    endpoints = endpoint_series(
        pd.read_csv(ENDPOINT, dtype={"NUTS_ID": str}))
    return annual_series(raw, endpoints, snow)


def validated_monthly(case: Setup):
    path = MONTHLY / f"{case.key}_2015_2024_monthly.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Expected municipal monthly observations at {path}. "
            "Run local_viirs_five_case_analysis.py first.")
    df = pd.read_csv(path)
    required = {"year", "calendar_month", "mean_radiance", "median_cf_cvg"}
    if required - set(df):
        raise ValueError(f"{path} missing {sorted(required - set(df))}")
    if len(df) != 120 or df.duplicated(["year", "calendar_month"]).any():
        raise ValueError(f"{path}: need exactly 120 unique 2015–24 months")
    if (set(zip(df.year, df.calendar_month)) !=
            {(year, mon) for year in range(2015, 2025)
             for mon in range(1, 13)}):
        raise ValueError(f"{path}: missing/out-of-window month")
    df["date"] = pd.to_datetime(dict(year=df.year, month=df.calendar_month,
                                     day=1))
    df = df.sort_values("date")
    df["mean_radiance"] = pd.to_numeric(df.mean_radiance, errors="raise")
    if not np.isfinite(df.mean_radiance).all():
        raise ValueError(f"{path}: nonfinite monthly radiance")
    return df


def baseline_and_anomaly(data: pd.DataFrame, case: Setup):
    baseline = (data.loc[data.year.isin(case.references)]
                .groupby("calendar_month").mean_radiance.mean())
    if len(baseline) != 12 or not np.isfinite(baseline).all():
        raise ValueError("Reference years must yield exactly 12 months")
    event = data.loc[data.year.eq(case.outcome_year)].copy()
    if len(event) != 12:
        raise ValueError("Event year must yield 12 months")
    event["reference_mean"] = event.calendar_month.map(baseline)
    event["anomaly"] = event.mean_radiance - event.reference_mean
    return event


def event_shading(ax, case: Setup, *, annotation=False):
    start, end = (pd.Timestamp(s) for s in case.event_window)
    ax.axvspan(start, end, color="#e1bb7c", alpha=.28, zorder=0)
    if case.led_window:
        a, b = (pd.Timestamp(s) for s in case.led_window)
        ax.axvspan(a, b, color="#90b7d3", alpha=.2, zorder=0)
    if annotation:
        # Identically dated rollout markers, not VIIRS-optimised changepoints.
        for i, (date, label, kind) in enumerate(case.dates):
            loc = pd.Timestamp(date)
            ax.axvline(loc, color="#875b32" if kind == "rollout" else "#57728e",
                       linestyle=":", linewidth=1, alpha=.75)
            if kind == "approval" or (case.key == "echirolles" and
                                      date in ("2022-07-29", "2022-10-01")):
                ax.annotate(label, xy=(loc, 1), xycoords=(
                    ax.get_xaxis_transform()),
                    xytext=(2, -9-18*(i % 2)), textcoords="offset points",
                    rotation=90, va="top", ha="left",
                    fontsize=7.7, color="#514b45")


def render(case: Setup, annual: pd.DataFrame, monthly: pd.DataFrame,
           destination: Path, *, formats=("png", "pdf", "svg")):
    r = annual.loc[annual.NUTS_ID.eq(case.region)].sort_values("year")
    if len(r) != 11 or list(r.year) != list(range(2014, 2025)):
        raise ValueError(f"Missing verified full annual NUTS3 for {case.region}")
    if not np.isfinite(r.radiance_locked_nw_cm2_sr).all():
        raise ValueError("Invalid annual radiance")
    event = baseline_and_anomaly(monthly, case)
    current = monthly[monthly.date >= pd.Timestamp(case.first_month)].copy()

    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.size": 9.5, "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlelocation": "left",
    })
    fig, axes = plt.subplots(
        3, 1, figsize=(12.8, 10.2), constrained_layout=True,
        gridspec_kw={"height_ratios": [1.0, 1.18, .85]})
    aa, am, ad = axes
    vals = r.radiance_locked_nw_cm2_sr.to_numpy(float)
    aa.axvspan(case.outcome_year-.45, case.outcome_year+.45,
               color="#e1bb7c", alpha=.25, zorder=0)
    if case.led_window:
        aa.axvspan(2023.48, 2024.48, color="#90b7d3", alpha=.18)
    aa.plot(r.year, vals, "-o", linewidth=1.9, markersize=4.8,
            color="#24465f", label="Locked annual VNL V2 NUTS2016 radiance")
    aa.set_xlim(2013.7, 2024.45)
    aa.set_xticks(list(range(2014, 2025)))
    aa.set_ylabel(f"Annual regional radiance\n({RUNIT})")
    aa.set_title(
        f"A   {case.region} | {REGIONS[case.region]['label']}, "
        "locked annual 2014–2024")
    aa.grid(axis="y", alpha=.2)
    aa.annotate(
        f"2014: {vals[0]:.2f}\n2024: {vals[-1]:.2f}"
        f"\nΔ {100*(vals[-1]/vals[0]-1):+.1f}%",
        xy=(.015, .06), xycoords="axes fraction", ha="left", va="bottom",
        fontsize=9, color="#24465f")
    aa.text(.99, .97, "2021: locked snow-corrected value",
            ha="right", va="top", transform=aa.transAxes,
            fontsize=8, color="#586675")
    aa.legend(loc="upper right", fontsize=8.5)

    am.plot(current.date, current.mean_radiance, "-",
            color="#969da2", linewidth=1.2, alpha=.8,
            label="Municipal monthly VNL V1, full period")
    selected = current.loc[current.year.eq(case.outcome_year)]
    am.plot(selected.date, selected.mean_radiance, "-o",
            markersize=4.2, linewidth=1.7, color="#193d5a",
            label=f"{case.outcome_year} municipality radiance", zorder=4)
    am.plot(event.date, event.reference_mean, "--",
            linewidth=1.3, color="#bf743c",
            label="Same-calendar-month mean, " + "/".join(
                map(str, case.references)))
    event_shading(am, case, annotation=True)
    am.set_xlim(pd.Timestamp(case.first_month),
                pd.Timestamp("2025-01-01"))
    am.set_ylabel(f"Municipal monthly radiance\n({RUNIT})")
    am.set_title(
        f"B   {case.title} municipality | monthly VNL V1 "
        f"({case.first_month[:4]}–2024)\n"
        f"       Tan: {case.event_label}"
        + ("  ·  Blue: subsequent LED programme" if case.led_window else ""))
    am.xaxis.set_major_locator(mdates.YearLocator())
    am.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    am.grid(axis="y", alpha=.22)
    am.legend(loc="upper left", ncol=2, fontsize=7.7, framealpha=.9)

    ad.plot(event.date, event.anomaly, "-o", color="#2c6080",
            linewidth=1.8, markersize=4)
    ad.axhline(0, color="#444", linewidth=.8)
    event_shading(ad, case)
    ad.set_xlim(pd.Timestamp(f"{case.outcome_year}-01-01"),
                pd.Timestamp(f"{case.outcome_year+1}-01-01"))
    ad.xaxis.set_major_locator(mdates.MonthLocator())
    ad.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ad.set_ylabel(f"Event-year monthly anomaly\n({RUNIT})")
    ad.set_title(
        f"C   {case.outcome_year} municipal monthly anomaly against "
        f"{'/'.join(map(str, case.references))} same-month reference")
    ad.grid(axis="y", alpha=.2)
    post = event[event.calendar_month.isin(case.comparison_months)]
    reference = post.reference_mean.mean()
    observed = post.mean_radiance.mean()
    ad.text(
        .99, .95,
        f"Selected post/late months: {','.join(map(str, case.comparison_months))}\n"
        f"Mean ΔR: {observed-reference:+.2f} ({100*(observed/reference-1):+.1f}%)",
        transform=ad.transAxes, ha="right", va="top",
        fontsize=9, bbox={"facecolor":"white", "alpha":.82,
                          "edgecolor":"#dddddd", "pad":4})

    fig.suptitle(
        f"{case.title}: regional history and independently dated municipal "
        "light operations\n"
        f"Annual: NUTS3 VNL V2  |  Monthly: municipality VNL V1  "
        "|  NOT numerically interchangeable",
        fontsize=13, weight="bold", y=1.045)
    fig.text(.02, -.035,
             case.notes + "  Annual observation year aggregates pre- "
             "and post-event months. No LED spectral correction or causal "
             "fraction estimated.",
             ha="left", va="top", fontsize=8.2, color="#555")
    destination.mkdir(parents=True, exist_ok=True)
    stem = f"{case.key}_nuts3_annual_and_municipal_monthly_2014_2024"
    for ext in formats:
        fig.savefig(destination / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return event, stem


def self_test():
    from export_coventry_echirolles_nuts3_bridge import self_test as bridge_test
    bridge_test()
    rng = np.random.default_rng(41)
    annual = pd.DataFrame([
        {"NUTS_ID": s.region, "year": year,
         "radiance_locked_nw_cm2_sr": (20 if s.region=="UKG33" else 3)
         - (year-2014)*(0.6 if s.region=="UKG33" else .1)}
        for s in CASES for year in range(2014, 2025)])
    months = pd.period_range("2015-01","2024-12",freq="M")
    with TemporaryDirectory(prefix="event-annual-monthly-") as tmp:
        for case in CASES:
            df = pd.DataFrame({
                "year": months.year, "calendar_month": months.month,
                "date": months.to_timestamp(),
                "mean_radiance": (20. + rng.normal(0, 1.4, len(months)))})
            e, stem = render(case, annual, df, Path(tmp), formats=("png",))
            assert len(e)==12 and (Path(tmp)/f"{stem}.png").is_file()
            assert all(np.isfinite(e.anomaly))
    print("Self-test PASS: read-only annual bridge + two full figure renders")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    annual = validated_annual()
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    annual[["NUTS_ID","region_name","year",
            "radiance_locked_nw_cm2_sr","snow_corrected_2021",
            "delta_previous_year_nw","delta_previous_year_pct"]].to_csv(
        DATA_OUT/"coventry_isere_plot_validated_annual.csv",index=False)
    for case in CASES:
        monthly = validated_monthly(case)
        event, stem = render(case,annual,monthly,OUT)
        event[["month","year","calendar_month","mean_radiance",
               "reference_mean","anomaly","median_cf_cvg"]].to_csv(
            DATA_OUT/f"{case.key}_event_year_monthly_anomalies.csv",
            index=False) if "month" in event else event[
                ["date","year","calendar_month","mean_radiance",
                 "reference_mean","anomaly","median_cf_cvg"]].to_csv(
            DATA_OUT/f"{case.key}_event_year_monthly_anomalies.csv",
            index=False)
        print(f"{case.title}: {OUT / (stem+'.png')}")
    print(f"Editable PDF/SVG plus PNG for both cases. Tables: {DATA_OUT}")
    print("All frozen parquet, source rasters, snow cache and PRC inputs READ ONLY.")


if __name__=="__main__":
    main()
