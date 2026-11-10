#!/usr/bin/env python3
"""Direct 2015→2019 measured sky BRIGHTNESS versus historical VIIRS.

Uses EXISTING output of plot_madrid_viirs_robles_2021.py: no raster re-read,
no recomputation/change of annual NUTS3 VIIRS, 2021 snow or PRC.

  python src/plot_madrid_brightness_viirs_direct_contrast.py --self-test
  python src/plot_madrid_brightness_viirs_direct_contrast.py

Bands are different photometric passbands; compare only WITHIN-band
2015-indexed changes. VIIRS ES300 ≠ VIIRS Madrid municipality ≠ UCM AstMon
zenith sky. Years 2011, 2015, 2019 are actual published Table 6 anchors,
not interpolated annual observations. No CI inferred from published spreads.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data/processed/local_lighting_interventions/madrid_ground_viirs"
FIG = ROOT / "figures/local_lighting_interventions/madrid_ground_viirs"
ANNUAL = DATA / "madrid_ES300_locked_annual_vnlv2_2014_2024.csv"
MONTHLY = DATA / "madrid_municipality_vcmslcfg_2015_2024.csv"
GROUND = DATA / "robles2021_table6_after_midnight_yearly_medians.csv"

COLORS = {
    "NUTS3 VIIRS": "#243f58",
    "City VIIRS": "#698c9e",
    "Johnson B": "#347cb7",
    "Johnson V": "#37895f",
    "Johnson R": "#b9564d",
}


def validate(annual, monthly, ground):
    if (annual.duplicated(["year"]).any() or len(annual) != 11 or
            set(annual.year.astype(int)) != set(range(2014, 2025))):
        raise ValueError("Need complete unique 2014–24 frozen ES300 annual CSV")
    if (monthly.duplicated(["year", "calendar_month"]).any() or
            len(monthly) != 120 or set(monthly.year.astype(int))
            != set(range(2015, 2025))):
        raise ValueError("Need 120 unique 2015–24 municipal monthly observations")
    for year in range(2015, 2025):
        if set(monthly.loc[monthly.year.eq(year), "calendar_month"]) != set(
                range(1, 13)):
            raise ValueError(f"Missing same-month observations {year}")
    if len(ground) != 15 or ground.duplicated(["year", "band"]).any():
        raise ValueError("Published Table 6 should have exactly 15 rows")
    if set(ground.year) != {2011, 2015, 2019} or set(ground.band) != {
            "B", "V", "R", "B-V", "V-R"}:
        raise ValueError("Unexpected Table 6 years or bands")
    for name, values in [
        ("annual", annual.radiance_locked),
        ("municipal", monthly.municipality_mean_vnlv1),
        ("ground", ground.published_median_mag_per_arcsec2),
    ]:
        if not np.isfinite(values.astype(float)).all():
            raise ValueError(f"Nonfinite {name}")
    if not (annual.radiance_locked.gt(0).all() and
            monthly.municipality_mean_vnlv1.gt(0).all()):
        raise ValueError("Unexpected non-positive VIIRS radiance")
    # Refuse column misinterpreted as flux where astronomical magnitudes expected
    expected = {2015: {"B": 19.27, "V": 18.07, "R": 17.42},
                2019: {"B": 18.98, "V": 17.95, "R": 17.44}}
    for year, bands in expected.items():
        for band, mag in bands.items():
            actual = ground.loc[
                ground.year.eq(year) & ground.band.eq(band),
                "published_median_mag_per_arcsec2"
            ].iloc[0]
            if not np.isclose(actual, mag, atol=.001):
                raise ValueError(f"Unexpected paper Table 6 {year} {band} magnitude")


def comparison(annual, monthly, ground):
    validate(annual, monthly, ground)
    regional = annual.set_index("year").radiance_locked
    local = monthly.groupby("year").municipality_mean_vnlv1.mean()
    rows = []
    for label, source, group in [
        ("NUTS3 VIIRS", regional, "satellite"),
        ("City VIIRS", local, "satellite"),
    ]:
        old, new = float(source.loc[2015]), float(source.loc[2019])
        rows.append({
            "measure": label, "measurement_type": group,
            "year_start": 2015, "year_end": 2019,
            "value_start": old, "value_end": new,
            "relative_change_pct": (new/old-1)*100,
            "method": "yearly annual VNL V2" if label == "NUTS3 VIIRS"
            else "mean of 12 VNL V1 monthly city means"})
    for band in ("B", "V", "R"):
        bb = ground.loc[ground.band.eq(band)].set_index("year")
        mag_start, mag_end = (float(bb.loc[year,
            "published_median_mag_per_arcsec2"]) for year in (2015, 2019))
        # The photometric zero point cancels in WITHIN-band ratio only.
        pct = (10 ** (-.4*(mag_end - mag_start)) - 1)*100
        rows.append({
            "measure": f"Johnson {band}", "measurement_type": "ground zenith",
            "year_start": 2015, "year_end": 2019,
            "value_start": mag_start, "value_end": mag_end,
            "relative_change_pct": pct,
            "method": "Robles et al. Table 6; 10^(-0.4*Δmag)"})
    return pd.DataFrame(rows)


def plot(annual, monthly, ground, folder, *, formats=("png", "pdf", "svg")):
    changes = comparison(annual, monthly, ground)
    regional = annual.set_index("year").radiance_locked.sort_index()
    local = monthly.groupby("year").municipality_mean_vnlv1.mean()
    plt.rcParams.update({
        "font.size": 9.4, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(
        3, 1, figsize=(12.8, 11),
        constrained_layout=True, gridspec_kw={
            "height_ratios": [1.05, .95, 1.1]})
    top, mid, bot = axes
    names = changes.measure.tolist()
    vals = changes.relative_change_pct.to_numpy(float)
    x = np.array([0, 1, 3, 4, 5])
    top.bar(x, vals, color=[COLORS[n] for n in names], width=.68)
    top.axhline(0, color="#424242", lw=.9)
    top.axvline(2, color="#888", lw=.9, ls=":", alpha=.7)
    top.set_xticks(x, [
        "ES300\nannual VIIRS", "Madrid\nmonthly VIIRS",
        "AstMon\nJohnson B", "AstMon\nJohnson V", "AstMon\nJohnson R"])
    top.set_ylabel("2015→2019 change (%)\npositive = brighter")
    top.set_title(
        "A  DIRECT RESULT: measured ground sky brightness can increase "
        "while VIIRS decreases", loc="left", weight="bold")
    top.grid(axis="y", alpha=.2)
    top.set_axisbelow(True)
    top.text(0.19, .95, "Upward satellite radiance\n(region / municipality)",
             transform=top.transAxes, ha="center", va="top",
             fontsize=8, color="#55616a")
    top.text(.78, .95, "Ground zenith brightness\n(one UCM observatory)",
             transform=top.transAxes, ha="center", va="top",
             fontsize=8, color="#55616a")
    for xpos, v in zip(x, vals):
        top.annotate(f"{v:+.1f}%", xy=(xpos, v),
                     xytext=(0, 5 if v >= 0 else -15),
                     textcoords="offset points", ha="center",
                     fontsize=9, weight="bold")

    mid.plot(regional.index, 100 * regional / regional.loc[2015],
             "-o", color=COLORS["NUTS3 VIIRS"],
             label="ES300 annual VNL V2", lw=1.8, ms=4.5)
    mid.plot(local.index, 100 * local / local.loc[2015],
             "-o", color=COLORS["City VIIRS"],
             label="City mean of monthly VNL V1", lw=1.65, ms=4.0)
    mid.set_xlim(2013.8, 2024.25)
    mid.set_xticks(np.arange(2014, 2025, 1))
    mid.axhline(100, color="#888", lw=.8, ls="--")
    mid.axvspan(2014.8, 2016.9, color="#dba568", alpha=.19)
    mid.set_ylabel("Relative VIIRS radiance\n2015 = 100 (within series)")
    mid.set_title("B  VIIRS historical radiance, two distinct geographies/products",
                  loc="left")
    mid.legend(loc="upper right", fontsize=8)
    mid.grid(axis="y", alpha=.2)

    offsets = {"B": -.09, "V": 0., "R": .09}
    for band, name in [("B", "Johnson B"), ("V", "Johnson V"),
                       ("R", "Johnson R")]:
        sub = ground.loc[ground.band.eq(band)].sort_values("year")
        base_mag = float(sub.loc[sub.year.eq(2015),
                                 "published_median_mag_per_arcsec2"].iloc[0])
        rel = 100 * 10**(-.4*(sub.published_median_mag_per_arcsec2-base_mag))
        bot.scatter(sub.year+offsets[band], rel,
                    color=COLORS[name], s=53, label=name,
                    edgecolor="white", linewidth=.45, zorder=3)
        for year, y in zip(sub.year, rel):
            if year in (2015, 2019):
                bot.annotate(f"{y:.1f}", (year+offsets[band], y),
                             xytext=(0, 6 if band!="R" else -12),
                             textcoords="offset points", ha="center",
                             fontsize=7, color=COLORS[name])
    bot.set_xticks([2011, 2015, 2019])
    bot.set_xlim(2010.45, 2019.55)
    bot.set_ylabel("Ground zenith sky brightness\n2015 = 100 within each band")
    bot.axhline(100, color="#888", lw=.8, ls="--")
    bot.axvspan(2014.8, 2016.9, color="#dba568", alpha=.19)
    bot.set_title(
        "C  UCM AstMon AFTER-midnight sky brightness: published 2011, "
        "2015, 2019 observations", loc="left")
    bot.legend(loc="upper center", ncol=3, fontsize=8)
    bot.grid(axis="y", alpha=.2)

    fig.suptitle(
        "Madrid: change in measured SKY BRIGHTNESS versus VIIRS radiance\n"
        "Ground zenith photometry and orbital upward radiance are not equivalent",
        fontsize=13.5, weight="bold", y=1.035)
    fig.text(
        .01, -.038,
        "Ground: Robles et al. (2021), Table 6; after-midnight, cloud/moon "
        "filtered; site at UCM. VIIRS: ES300 frozen annual 2014–24 and Madrid "
        "municipality monthly means. Shades: broad 2015–16 retrofit. "
        "2011/2015/2019 are discrete measured ground anchors, NOT interpolated. "
        "No causal LED/optical/dimming partition inferred.",
        fontsize=7.6, color="#555", va="top")
    folder.mkdir(parents=True, exist_ok=True)
    name = "madrid_2015_2019_direct_sky_brightness_vs_viirs"
    for ext in formats:
        fig.savefig(folder/f"{name}.{ext}", dpi=190, bbox_inches="tight")
    plt.close(fig)
    return changes, name


def self_test():
    annual = pd.DataFrame({
        "year": range(2014, 2025),
        "radiance_locked": np.linspace(12., 8., 11)})
    m = pd.period_range("2015-01", "2024-12", freq="M")
    monthly = pd.DataFrame({
        "year": m.year, "calendar_month": m.month,
        "municipality_mean_vnlv1": 30 + np.sin(np.arange(120)/6)})
    table = {
        2011: {"B":19.15, "V":17.82, "R":17.11,
               "B-V":1.32,"V-R":.78},
        2015: {"B":19.27, "V":18.07, "R":17.42,
               "B-V":1.20,"V-R":.66},
        2019: {"B":18.98, "V":17.95, "R":17.44,
               "B-V":1.02,"V-R":.51}}
    ground = pd.DataFrame([
        {"year":year,"band":band,
         "published_median_mag_per_arcsec2":val}
        for year,bands in table.items() for band,val in bands.items()])
    with TemporaryDirectory(prefix="madrid-brightness-") as folder:
        changes, name = plot(annual, monthly, ground, Path(folder))
        assert all((Path(folder)/f"{name}.{ext}").exists()
                   for ext in ("png","pdf","svg"))
        b = float(changes.loc[changes.measure.eq("Johnson B"),
                              "relative_change_pct"].iloc[0])
        assert np.isclose(b, 30.62, atol=.2), b
    print("Self-test PASS: 2015→19 sky-brightness conversion, source "
          "identity guards, PNG/PDF/SVG figure rendering")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    missing = [p for p in (ANNUAL, MONTHLY, GROUND) if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "First run python src/plot_madrid_viirs_robles_2021.py "
            "to generate these verified inputs: "
            + ", ".join(str(p) for p in missing))
    a, m, g = pd.read_csv(ANNUAL), pd.read_csv(MONTHLY), pd.read_csv(GROUND)
    changes, name = plot(a, m, g, FIG)
    changes.to_csv(DATA/"madrid_2015_2019_ground_vs_viirs_brightness.csv",
                   index=False)
    print(changes.to_string(index=False))
    print(f"Figure: {FIG / (name+'.png')}")
    print("Read-only derived figure. Source monthly rasters, annual parquet, "
          "2021 snow patch and PRC are unchanged.")


if __name__ == "__main__":
    main()
