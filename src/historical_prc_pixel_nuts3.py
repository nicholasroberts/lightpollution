#!/usr/bin/env python3
"""Historical PRC from individual EOG VIIRS pixels, aggregated AFTER transformation.

Annual and 2014-relative changes in area-weighted NUTS2016 means; maps and
NUTS3 rankings. Native-resolution pixel rasters may also be written for
selected years. The NUTS3 2021 snow adjustment is NOT a pixel correction.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window, from_bounds, transform as window_transform
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
import yaml

from extract_nuts3_radiance import discover_rasters, nuts_path_for_release, validate_raster
from project_config import PROJECT_ROOT, load_config


def prc_of_radiance(raw, offset, yf, y0, k):
    """Transform masked VIIRS values; zero is valid, negative radiance invalid."""
    a = np.ma.asarray(raw)
    values = np.asarray(a.data, dtype=np.float32)
    invalid = np.ma.getmaskarray(a) | ~np.isfinite(values) | (values < 0)
    r = np.maximum(values - offset, 0.0)
    result = yf + (y0 - yf) * np.exp(-k * r)
    return np.ma.masked_array(result.astype(np.float32), mask=invalid)


class PRCRasterSource(RasterioRasterSource):
    def __init__(self, ds, offset, params):
        super().__init__(ds, 1, name="prc")
        self.offset = offset
        self.params = params

    def read_window(self, x0, y0, nx, ny):
        return prc_of_radiance(super().read_window(x0, y0, nx, ny),
                               self.offset, *self.params)


def load_inputs(config, transfer_path):
    with transfer_path.open(encoding="utf-8") as f:
        transfer = yaml.safe_load(f)
    assert transfer["model_family"] == "R_SSasymp_exponential"
    params = transfer["parameters"]
    p = (float(params["yf"]), float(params["y0"]), float(params["k"]))
    offset = float(config["calibration"]["pixel_offset_nW_cm2_sr"])
    from_year = int(config["calibration"]["apply_from_year"])
    expected = transfer["radiance"]["pixel_calibration"]
    if (abs(offset - float(expected["positive_pixel_subtraction"])) > 1e-9 or
            from_year != int(expected["apply_from_year"])):
        raise ValueError("PRC config calibration does not match pipeline calibration")
    if not (0 <= p[0] < p[1] and p[2] > 0):
        raise ValueError("Unexpected PRC function parameters")
    return transfer, p, offset, from_year


def load_nuts(config):
    path = nuts_path_for_release(config, 2016)
    if not path.is_file():
        raise FileNotFoundError(path)
    nuts = gpd.read_file(path)
    if nuts.crs is None:
        raise ValueError("NUTS polygons lack CRS")
    nuts = nuts.to_crs("EPSG:4326")
    if nuts["NUTS_ID"].duplicated().any():
        raise ValueError("Duplicate NUTS2016 ids")
    return nuts


def extract_year(path, year, nuts, config, params, offset, from_year):
    validate_raster(path, config["nuts"]["crs"])
    off = offset if year >= from_year else 0.0
    print(f"PRC for {year}: {path.name}; VIIRS pixel correction={off:g}",
          flush=True)
    includes = [v for v in ["NUTS_ID", "CNTR_CODE", "NUTS_NAME", "NAME_LATN"]
                if v in nuts.columns]
    with rasterio.open(path) as ds:
        source = PRCRasterSource(ds, off, params)
        extracted = exact_extract(
            source, nuts,
            ["mean(coverage_weight=area_spherical_m2)",
             "count(coverage_weight=area_spherical_km2)"],
            include_cols=includes, output="pandas", progress=True)
    # exactextract labels single input using "prc_" when a named RasterSource is used
    rename = {}
    for target, options in [
        ("prc_mean", ["prc_mean", "mean"]),
        ("valid_area_km2", ["prc_count", "count"]),
    ]:
        found = [v for v in options if v in extracted.columns]
        if len(found) != 1:
            raise ValueError(f"Cannot identify {target} in {list(extracted.columns)}")
        rename[found[0]] = target
    extracted = extracted.rename(columns=rename)
    extracted["year"] = year
    extracted["viirs_file"] = path.name
    extracted["pixel_offset"] = off
    return extracted


def build_panel(d, years, min_coverage):
    if d.duplicated(["NUTS_ID", "year"]).any():
        raise ValueError("Duplicate region-years")
    d = d.sort_values(["NUTS_ID", "year"]).reset_index(drop=True)
    base = d.loc[d.year.eq(2014)].set_index("NUTS_ID")
    if len(base) == 0:
        raise ValueError("2014 baseline missing")
    d["prc_2014"] = d.NUTS_ID.map(base.prc_mean)
    d["valid_area_2014_km2"] = d.NUTS_ID.map(base.valid_area_km2)
    d["coverage_ratio_to_2014"] = d.valid_area_km2 / d.valid_area_2014_km2
    grp = d.groupby("NUTS_ID", sort=False)
    d["previous_year"] = grp.year.shift()
    d["previous_prc"] = grp.prc_mean.shift()
    d["annual_prc_change_absolute"] = (d.prc_mean - d.previous_prc).where(
        d.year.sub(d.previous_year).eq(1))
    d["annual_prc_change_pct"] = (
        100.0 * d.annual_prc_change_absolute / d.previous_prc)
    d["cumulative_prc_change_absolute"] = d.prc_mean - d.prc_2014
    d["cumulative_prc_change_pct"] = (
        100.0 * d.cumulative_prc_change_absolute / d.prc_2014)
    # A gross coverage change can make different-year polygon means misleading.
    # Require comparable valid areas in EVERY year of the requested interval
    # for headline rankings; this does not establish exact common-pixel support.
    g = d.groupby("NUTS_ID")
    d["coverage_stable"] = (
        d.NUTS_ID.map(g.coverage_ratio_to_2014.min()) >= min_coverage
    ) & (d.NUTS_ID.map(g.coverage_ratio_to_2014.max()) <= 1 / min_coverage)
    d["coverage_stable"] &= (d.valid_area_2014_km2 > 0)
    d["coverage_stable"] &= d.NUTS_ID.map(g.year.nunique()).eq(len(years))
    d["coverage_stable"] &= np.isfinite(d.prc_mean) & np.isfinite(d.prc_2014)
    return d


def europe_table(d):
    rows = []
    for year, x in d.groupby("year"):
        x = x[np.isfinite(x.prc_mean) & np.isfinite(x.valid_area_km2)
              & (x.valid_area_km2 > 0)]
        rows.append({"year": int(year), "nuts3_regions": len(x),
                     "nuts_covered_valid_km2": x.valid_area_km2.sum(),
                     "area_weighted_prc": np.average(x.prc_mean,
                                                       weights=x.valid_area_km2)})
    out = pd.DataFrame(rows).sort_values("year")
    out["annual_prc_change_absolute"] = out.area_weighted_prc.diff()
    out["annual_prc_change_pct"] = out.area_weighted_prc.pct_change() * 100
    out["cumulative_prc_change_absolute"] = (
        out.area_weighted_prc - out.loc[out.year.eq(2014), "area_weighted_prc"].iloc[0])
    out["cumulative_prc_change_pct"] = (
        100 * out.cumulative_prc_change_absolute /
        out.loc[out.year.eq(2014), "area_weighted_prc"].iloc[0])
    return out


METRICS = [
    ("annual_prc_change_absolute", "Annual absolute PRC change", "annual_absolute_prc"),
    ("annual_prc_change_pct", "Annual PRC change (%)", "annual_percent_prc"),
    ("cumulative_prc_change_absolute", "2014-relative absolute PRC change",
     "cumulative_absolute_prc"),
    ("cumulative_prc_change_pct", "2014-relative PRC change (%)",
     "cumulative_percent_prc"),
]


def plot_region_atlases(d, nuts, output, years):
    mainland = nuts.copy()
    common = d[[ "NUTS_ID", "year", *[m[0] for m in METRICS]]]
    bounds = (-12, 36, 34, 73)
    for metric, label, stem in METRICS:
        all_values = common[metric].to_numpy()
        finite = np.abs(all_values[np.isfinite(all_values)])
        vmax = max(float(np.quantile(finite, .98)) if finite.size else 0, 1e-4)
        with PdfPages(output / f"{stem}_nuts3_atlas.pdf") as pdf:
            for year in years:
                if year == 2014:
                    continue
                this = common.loc[common.year.eq(year), ["NUTS_ID", metric]]
                g = mainland.merge(this, on="NUTS_ID", how="left")
                fig, ax = plt.subplots(figsize=(10.5, 8.5))
                g.plot(column=metric, ax=ax, cmap="RdBu", vmin=-vmax, vmax=vmax,
                       linewidth=0, missing_kwds={"color": "lightgray"},
                       legend=True, legend_kwds={"label": label, "shrink": .65})
                ax.set(xlim=(bounds[0], bounds[2]), ylim=(bounds[1], bounds[3]),
                       title=f"{label}: {year}\n"
                             f"Pixel-level PRC transformed before NUTS2016 averaging")
                ax.set_xlabel("Longitude")
                ax.set_ylabel("Latitude")
                ax.text(.02, .02,
                        "Display: geographic mainland window; tables retain all NUTS2016 regions.",
                        transform=ax.transAxes, fontsize=8,
                        bbox={"facecolor": "white", "alpha": .85, "edgecolor": "none"})
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
        print("Wrote", output / f"{stem}_nuts3_atlas.pdf", flush=True)


def plot_europe(e, output):
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), layout="constrained")
    for ax, (col, name, _) in zip(axs.flat, METRICS):
        ax.plot(e.year, e[col], marker="o")
        ax.axhline(0, color="gray", linestyle="--", linewidth=.75)
        ax.set(title=name, xlabel="Year", ylabel="PRC units" if "absolute" in col else "%")
        ax.grid(alpha=.2)
        ax.set_xticks(e.year)
        ax.tick_params(axis="x", labelrotation=45)
    fig.suptitle("Area-weighted PRC changes over the NUTS2016-covered area")
    fig.savefig(output / "prc_europe_annual_cumulative_changes.pdf")
    plt.close(fig)


def write_pixel_year(year, rasters, params, offset, from_year, output,
                     extent=(-12, 34, 36, 73)):
    """Optional native-resolution 5-band VIIRS pixel PRC and change GeoTIFF."""
    if year <= 2014 or year - 1 not in rasters:
        raise ValueError("Pixel output year requires 2014 and immediately preceding year")
    paths = [rasters[2014], rasters[year - 1], rasters[year]]
    with ExitStack() as stack:
        ds = [stack.enter_context(rasterio.open(p)) for p in paths]
        if (any(x.crs != ds[0].crs or x.width != ds[0].width or
                x.height != ds[0].height or not x.transform.almost_equals(ds[0].transform)
                for x in ds[1:])):
            raise ValueError("Original pixel rasters do not share grid; regridding required")
        if ds[0].crs.to_epsg() != 4326:
            raise ValueError("Pixel writer requires EPSG:4326")
        full = from_bounds(*extent, transform=ds[0].transform)
        full = full.round_offsets().round_lengths()
        if full.col_off < 0 or full.row_off < 0 or (
            full.col_off + full.width > ds[0].width or
            full.row_off + full.height > ds[0].height):
            raise ValueError("Pixel export window outside original VIIRS grid")
        profile = ds[0].profile.copy()
        profile.update(driver="GTiff", height=int(full.height), width=int(full.width),
                       count=5, dtype="float32",
                       transform=window_transform(full, ds[0].transform),
                       nodata=np.nan, tiled=True, blockxsize=512, blockysize=512,
                       compress="deflate", predictor=3, BIGTIFF="IF_SAFER")
        path = output / f"prc_pixels_{year}_annual_and_2014_cumulative.tif"
        with rasterio.open(path, "w", **profile) as dst:
            for idx, desc in enumerate([
                f"PRC_{year}", "annual_delta_prc", "annual_delta_pct",
                "from_2014_delta_prc", "from_2014_delta_pct"], start=1):
                dst.set_band_description(idx, desc)
            for start in range(0, int(full.height), 512):
                h = min(512, int(full.height) - start)
                src_window = Window(full.col_off, full.row_off + start,
                                    full.width, h)
                local = Window(0, start, full.width, h)
                vals = [prc_of_radiance(
                    x.read(1, window=src_window, masked=True),
                    offset if yr >= from_year else 0.0, *params)
                    for x, yr in zip(ds, [2014, year - 1, year])]
                base, prev, cur = vals
                paired_annual = np.ma.getmaskarray(prev) | np.ma.getmaskarray(cur)
                paired_cum = np.ma.getmaskarray(base) | np.ma.getmaskarray(cur)
                def band(x, invalid):
                    return np.ma.masked_array(x, mask=invalid).filled(np.nan).astype("float32")
                arrays = [
                    band(cur, np.ma.getmaskarray(cur)),
                    band(cur - prev, paired_annual),
                    band(100.0 * (cur - prev) / prev, paired_annual),
                    band(cur - base, paired_cum),
                    band(100.0 * (cur - base) / base, paired_cum),
                ]
                for b, values in enumerate(arrays, start=1):
                    dst.write(values, b, window=local)
                if start % 4096 == 0:
                    print(f"  Native pixels {year}: {start}/{full.height} rows", flush=True)
        print("Wrote native-resolution raster:", path, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--transfer-config", default="config/prc_transfer.yaml")
    parser.add_argument("--pilot", action="store_true",
                        help="2014/2018/2019/2021/2024 preliminary extraction")
    parser.add_argument("--pixel-year", action="append", type=int, default=[],
                        help="Optional native raster (e.g. --pixel-year 2024)")
    parser.add_argument("--no-maps", action="store_true")
    parser.add_argument("--min-coverage", type=float, default=.95)
    args = parser.parse_args()
    if not (0.5 <= args.min_coverage <= 1):
        raise ValueError("--min-coverage must be between .5 and 1")
    config = load_config(args.config)
    transfer_path = Path(args.transfer_config)
    if not transfer_path.is_absolute():
        transfer_path = PROJECT_ROOT / transfer_path
    transfer, params, offset, from_year = load_inputs(config, transfer_path)
    years = [2014, 2018, 2019, 2021, 2024] if args.pilot else list(range(2014, 2025))
    rasters = discover_rasters(config, years)
    nuts = load_nuts(config)
    output = PROJECT_ROOT / "data/processed/historical_prc_2014_2024"
    figures = PROJECT_ROOT / "figures/historical_prc_2014_2024"
    output.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    parts = [extract_year(rasters[y], y, nuts, config, params, offset, from_year)
             for y in years]
    panel = build_panel(pd.concat(parts, ignore_index=True),
                        years, args.min_coverage)
    prefix = "pilot_" if args.pilot else ""
    panel.to_csv(output / (prefix + "prc_nuts3_annual_cumulative.csv"), index=False)
    panel.to_parquet(output / (prefix + "prc_nuts3_annual_cumulative.parquet"),
                     index=False)
    europe = europe_table(panel)
    europe.to_csv(output / (prefix + "prc_europe_annual_cumulative.csv"), index=False)
    last = panel.loc[panel.year.eq(2024)].copy()
    last = last.loc[last.coverage_stable &
                    np.isfinite(last.cumulative_prc_change_absolute)]
    last["rank_largest_gain"] = last.cumulative_prc_change_absolute.rank(
        ascending=False, method="min")
    last["rank_largest_loss"] = last.cumulative_prc_change_absolute.rank(
        ascending=True, method="min")
    last["rank_pct_gain"] = last.cumulative_prc_change_pct.rank(
        ascending=False, method="min")
    last["rank_pct_loss"] = last.cumulative_prc_change_pct.rank(
        ascending=True, method="min")
    last.sort_values("cumulative_prc_change_absolute",
                     ascending=False).to_csv(
        output / (prefix + "prc_nuts3_rankings_2014_2024.csv"), index=False)
    if not args.no_maps:
        plot_region_atlases(panel, nuts, figures, years)
        plot_europe(europe, figures)
    for y in args.pixel_year:
        if args.pilot:
            raise ValueError("--pixel-year cannot be used with --pilot")
        write_pixel_year(y, rasters, params, offset, from_year, output)
    print("\nOutputs:", output, "and", figures)
    print("NUTS3 regions/year:", panel.groupby("year").NUTS_ID.nunique().to_dict())
    print("Ranked 2014->2024 coverage-stable regions:", len(last))
    print("Top PRC gains:")
    print(last.nlargest(10, "cumulative_prc_change_absolute")[
        ["NUTS_ID", "NUTS_NAME", "cumulative_prc_change_absolute",
         "cumulative_prc_change_pct"]].to_string(index=False))
    print("Top PRC losses:")
    print(last.nsmallest(10, "cumulative_prc_change_absolute")[
        ["NUTS_ID", "NUTS_NAME", "cumulative_prc_change_absolute",
         "cumulative_prc_change_pct"]].to_string(index=False))
    print("PRC function:", transfer["equation"], params)
    print("2021 is NOT corrected for the known pixel-scale snow/albedo anomaly.")


if __name__ == "__main__":
    main()
