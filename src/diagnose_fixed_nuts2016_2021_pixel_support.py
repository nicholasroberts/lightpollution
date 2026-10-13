#!/usr/bin/env python3
"""Test whether the 2021 VIIRS mask expansion is transient at pixel level.

This is a targeted mechanistic diagnostic using only the existing 2020-2023
raw EOG average_masked rasters and the locked fixed-NUTS2016 cohort.

For 2020->2021 every pixel is classified as:
- newly positive: 2020 <= 0 and 2021 > 0;
- lost positive: 2020 > 0 and 2021 <= 0;
- persistent positive: 2020 > 0 and 2021 > 0.

The regional 2021 radiance change is then decomposed exactly as:
    new-positive contribution
  + lost-positive contribution
  + persistent-positive radiance change.

Newly positive 2021 pixels are additionally followed into 2022 to determine
whether they immediately disappear again. Radiance change is also measured on
pixels that remain positive in all four years 2020-2023, providing a clean
common-support test of the 2021 intensity pulse and its 2022 reversal.

No production files are modified.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from exactextract import exact_extract
from exactextract.raster import RasterioRasterSource
from shapely.geometry import box

from diagnose_fixed_nuts2016_2021_reversal import load_locked_panel
from diagnose_fixed_nuts2016_mask_occupancy import attach_fixed_brightness
from extract_nuts3_radiance import discover_rasters, validate_raster
from project_config import configured_path, load_config
from validate_fixed_nuts2016_2014_2024 import load_fixed_nuts


YEARS = [2020, 2021, 2022, 2023]
MODES = [
    "new_2021_contribution",
    "lost_2021_contribution",
    "new_2021_area",
    "transient_new_2021_area",
    "transient_new_2021_contribution",
    "stable4_area",
    "stable4_change_2021",
    "stable4_change_2022",
]
Q_LABELS = ["Q1 darkest", "Q2", "Q3", "Q4", "Q5 brightest"]


def _as_data_valid(arr):
    if np.ma.isMaskedArray(arr):
        data = np.asarray(arr.data, dtype=np.float32)
        valid = ~np.ma.getmaskarray(arr)
    else:
        data = np.asarray(arr, dtype=np.float32)
        valid = np.ones(data.shape, dtype=bool)
    valid &= np.isfinite(data)
    return data, valid


class FourYearStateSource(RasterioRasterSource):
    """Derived pixel state/value source sharing four aligned annual rasters."""

    def __init__(self, datasets: dict[int, object], *, mode: str, name: str):
        if mode not in MODES:
            raise ValueError(f"Unknown mode: {mode}")
        super().__init__(datasets[2021], 1, name=name)
        self.mode = mode
        self.sources = {
            year: RasterioRasterSource(ds, 1, name=f"y{year}")
            for year, ds in datasets.items()
        }

    def read_window(self, x0, y0, nx, ny):
        arrays = {}
        valid_all = None
        for year, src in self.sources.items():
            data, valid = _as_data_valid(src.read_window(x0, y0, nx, ny))
            arrays[year] = data
            valid_all = valid if valid_all is None else (valid_all & valid)

        v20, v21, v22, v23 = (arrays[y] for y in YEARS)
        p20, p21, p22, p23 = (arrays[y] > 0 for y in YEARS)

        new21 = (~p20) & p21
        lost21 = p20 & (~p21)
        transient_new21 = new21 & (~p22)
        stable4 = p20 & p21 & p22 & p23

        out = np.zeros(v21.shape, dtype=np.float32)

        if self.mode == "new_2021_contribution":
            out[new21] = v21[new21]
        elif self.mode == "lost_2021_contribution":
            out[lost21] = -v20[lost21]
        elif self.mode == "new_2021_area":
            out[new21] = 1.0
        elif self.mode == "transient_new_2021_area":
            out[transient_new21] = 1.0
        elif self.mode == "transient_new_2021_contribution":
            out[transient_new21] = v21[transient_new21]
        elif self.mode == "stable4_area":
            out[stable4] = 1.0
        elif self.mode == "stable4_change_2021":
            out[stable4] = v21[stable4] - v20[stable4]
        elif self.mode == "stable4_change_2022":
            out[stable4] = v22[stable4] - v21[stable4]
        else:  # pragma: no cover
            raise RuntimeError(self.mode)

        if valid_all is not None and not valid_all.all():
            return np.ma.masked_array(out, mask=~valid_all, copy=False)
        return out


def run_exact_extract(datasets, nuts, id_field, *, progress):
    sources = [
        FourYearStateSource(datasets, mode=mode, name=mode)
        for mode in MODES
    ]
    result = exact_extract(
        sources,
        nuts,
        ["mean(coverage_weight=area_spherical_m2)"],
        include_cols=[id_field],
        output="pandas",
        progress=progress,
    )
    expected = {id_field} | {f"{mode}_mean" for mode in MODES}
    missing = expected - set(result.columns)
    if missing:
        raise RuntimeError(
            "Unexpected exactextract schema; missing "
            f"{sorted(missing)}; found {list(result.columns)}"
        )
    result = result.rename(
        columns={f"{mode}_mean": mode for mode in MODES}
    )
    return result[[id_field] + MODES].copy()


def integration_self_test():
    """Cheap real exactextract test of the custom four-year raster source."""
    arrays = {
        2020: np.array([[0.0, 1.0], [2.0, 3.0]], dtype=np.float32),
        2021: np.array([[4.0, 0.0], [3.0, 4.0]], dtype=np.float32),
        2022: np.array([[0.0, 0.0], [3.5, 4.5]], dtype=np.float32),
        2023: np.array([[0.0, 0.0], [4.0, 5.0]], dtype=np.float32),
    }
    memories = {}
    datasets = {}
    try:
        for year, arr in arrays.items():
            mem = MemoryFile()
            ds = mem.open(
                driver="GTiff",
                height=2,
                width=2,
                count=1,
                dtype="float32",
                crs="EPSG:4326",
                transform=from_origin(0, 2, 1, 1),
                nodata=-9999.0,
            )
            ds.write(arr, 1)
            memories[year] = mem
            datasets[year] = ds

        gdf = gpd.GeoDataFrame(
            {"NUTS_ID": ["TEST"]},
            geometry=[box(0, 0, 2, 2)],
            crs="EPSG:4326",
        )
        out = run_exact_extract(
            datasets, gdf, "NUTS_ID", progress=False
        ).iloc[0]

        if not (out["new_2021_contribution"] > 0):
            raise RuntimeError("Self-test failed: new-positive contribution.")
        if not (out["lost_2021_contribution"] < 0):
            raise RuntimeError("Self-test failed: lost-positive contribution.")
        if not np.isclose(
            out["transient_new_2021_area"],
            out["new_2021_area"],
            rtol=1e-6,
            atol=1e-8,
        ):
            raise RuntimeError("Self-test failed: transient 2021 entrant area.")
        if not np.isclose(
            out["transient_new_2021_contribution"],
            out["new_2021_contribution"],
            rtol=1e-6,
            atol=1e-8,
        ):
            raise RuntimeError(
                "Self-test failed: transient 2021 entrant contribution."
            )
        if not (out["stable4_area"] > 0):
            raise RuntimeError("Self-test failed: stable-four-year area.")
        if not (out["stable4_change_2021"] > 0):
            raise RuntimeError("Self-test failed: stable 2021 intensity change.")
        if not (out["stable4_change_2022"] > 0):
            raise RuntimeError("Self-test failed: stable 2022 intensity change.")
    finally:
        for ds in datasets.values():
            ds.close()
        for mem in memories.values():
            mem.close()

    print("2021 pixel-support integration self-test: PASS", flush=True)


def extraction_manifest(rasters: dict[int, Path]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "year": year,
                "viirs_file": rasters[year].name,
                "size_bytes": int(rasters[year].stat().st_size),
            }
            for year in YEARS
        ]
    )


def extract_or_load(panel: pd.DataFrame, config: dict):
    processed = configured_path(config, "processed")
    cache = processed / "fixed_nuts2016_2021_pixel_support_cache.parquet"
    manifest_path = (
        processed / "fixed_nuts2016_2021_pixel_support_cache_manifest.csv"
    )

    rasters = discover_rasters(config, YEARS)
    manifest = extraction_manifest(rasters)

    if cache.exists() and manifest_path.exists():
        cached_manifest = pd.read_csv(manifest_path)
        if cached_manifest.equals(manifest):
            print(f"Using cached 2021 pixel-support extraction: {cache}", flush=True)
            out = pd.read_parquet(cache)
            return out, manifest

    release = int(config["analysis"]["fixed_nuts2016_2014_2024"]["nuts_release"])
    nuts = load_fixed_nuts(config, release)
    keep_ids = set(panel["NUTS_ID"].astype(str))
    id_field = config["nuts"]["id_field"]

    for path in rasters.values():
        validate_raster(path, config["nuts"]["crs"])

    print(
        "Extracting 2020-2023 pixel-support states on fixed NUTS2016 geography",
        flush=True,
    )
    datasets = {}
    try:
        for year in YEARS:
            datasets[year] = rasterio.open(rasters[year])
        out = run_exact_extract(
            datasets, nuts, id_field, progress=True
        )
    finally:
        for ds in datasets.values():
            ds.close()

    out[id_field] = out[id_field].astype(str)
    out = out[out[id_field].isin(keep_ids)].copy()
    if out[id_field].duplicated().any():
        raise RuntimeError("Duplicate NUTS_ID in pixel-support extraction.")
    if set(out[id_field]) != keep_ids:
        raise RuntimeError("Pixel-support extraction does not match locked cohort.")

    out.to_parquet(cache, index=False)
    manifest.to_csv(manifest_path, index=False)
    print(f"Cached 2021 pixel-support extraction: {cache}", flush=True)
    return out, manifest


def build_regional(panel: pd.DataFrame, components: pd.DataFrame) -> pd.DataFrame:
    p20 = panel[panel["year"] == 2020].copy()
    p21 = panel[panel["year"] == 2021].copy()

    cols = [
        "NUTS_ID",
        "CNTR_CODE",
        "density_band",
        "brightness_quintile",
        "baseline_raw_nw",
    ]
    meta = p20[cols].copy()
    meta["raw_2020_nw"] = p20["radiance_mean_raw"].to_numpy(float)
    meta = meta.merge(
        p21[["NUTS_ID", "radiance_mean_raw"]].rename(
            columns={"radiance_mean_raw": "raw_2021_nw"}
        ),
        on="NUTS_ID",
        how="left",
        validate="one_to_one",
    )
    d = meta.merge(
        components,
        on="NUTS_ID",
        how="left",
        validate="one_to_one",
    )

    d["total_2021_change_nw"] = d["raw_2021_nw"] - d["raw_2020_nw"]
    d["persistent_2021_change_nw"] = (
        d["total_2021_change_nw"]
        - d["new_2021_contribution"]
        - d["lost_2021_contribution"]
    )
    d["decomposition_residual_nw"] = (
        d["total_2021_change_nw"]
        - d["new_2021_contribution"]
        - d["lost_2021_contribution"]
        - d["persistent_2021_change_nw"]
    )

    d["surviving_new_2021_area"] = (
        d["new_2021_area"] - d["transient_new_2021_area"]
    )
    d["surviving_new_2021_contribution"] = (
        d["new_2021_contribution"]
        - d["transient_new_2021_contribution"]
    )

    d["transient_fraction_of_new_area"] = np.where(
        d["new_2021_area"] > 0,
        d["transient_new_2021_area"] / d["new_2021_area"],
        np.nan,
    )
    d["transient_fraction_of_new_contribution"] = np.where(
        d["new_2021_contribution"] > 0,
        d["transient_new_2021_contribution"]
        / d["new_2021_contribution"],
        np.nan,
    )
    d["stable4_conditional_change_2021_nw"] = np.where(
        d["stable4_area"] > 0,
        d["stable4_change_2021"] / d["stable4_area"],
        np.nan,
    )
    d["stable4_conditional_change_2022_nw"] = np.where(
        d["stable4_area"] > 0,
        d["stable4_change_2022"] / d["stable4_area"],
        np.nan,
    )

    # Compare the direct pixel-state decomposition to the locked zonal means.
    # The persistent term is the exact remainder; new/lost contributions must
    # nevertheless be finite and have the correct signs.
    finite_cols = [
        "new_2021_contribution",
        "lost_2021_contribution",
        "new_2021_area",
        "transient_new_2021_area",
        "stable4_area",
        "stable4_change_2021",
        "stable4_change_2022",
    ]
    if (~np.isfinite(d[finite_cols].to_numpy(float))).any():
        raise RuntimeError("Non-finite pixel-support component.")
    if (d["new_2021_contribution"] < -1e-8).any():
        raise RuntimeError("New-positive contribution became negative.")
    if (d["lost_2021_contribution"] > 1e-8).any():
        raise RuntimeError("Lost-positive contribution became positive.")
    if (
        d["transient_new_2021_area"]
        > d["new_2021_area"] + 1e-8
    ).any():
        raise RuntimeError("Transient entrant area exceeds new entrant area.")

    return d


def summarize_group(d: pd.DataFrame, group: str) -> pd.DataFrame:
    rows = []
    for label, g in d.groupby(group, observed=True):
        new_area = float(g["new_2021_area"].median())
        transient_area = float(g["transient_new_2021_area"].median())
        new_contrib = float(g["new_2021_contribution"].median())
        transient_contrib = float(
            g["transient_new_2021_contribution"].median()
        )
        rows.append(
            {
                group: str(label),
                "n_nuts3": int(len(g)),
                "median_total_2021_change_nw": float(
                    g["total_2021_change_nw"].median()
                ),
                "median_new_positive_contribution_nw": new_contrib,
                "median_lost_positive_contribution_nw": float(
                    g["lost_2021_contribution"].median()
                ),
                "median_persistent_positive_change_nw": float(
                    g["persistent_2021_change_nw"].median()
                ),
                "median_new_positive_area_pp": 100.0 * new_area,
                "median_transient_new_area_pp": 100.0 * transient_area,
                "median_transient_fraction_of_new_area": float(
                    g["transient_fraction_of_new_area"].median()
                ),
                "median_transient_fraction_of_new_contribution": float(
                    g["transient_fraction_of_new_contribution"].median()
                ),
                "median_stable4_area_pct": float(
                    100.0 * g["stable4_area"].median()
                ),
                "median_stable4_conditional_change_2021_nw": float(
                    g["stable4_conditional_change_2021_nw"].median()
                ),
                "median_stable4_conditional_change_2022_nw": float(
                    g["stable4_conditional_change_2022_nw"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def country_summary(d: pd.DataFrame) -> pd.DataFrame:
    return summarize_group(d, "CNTR_CODE")


def overall_summary(d: pd.DataFrame, country: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "median_total_2021_change_nw",
        "median_new_positive_contribution_nw",
        "median_lost_positive_contribution_nw",
        "median_persistent_positive_change_nw",
        "median_new_positive_area_pp",
        "median_transient_new_area_pp",
        "median_transient_fraction_of_new_area",
        "median_transient_fraction_of_new_contribution",
        "median_stable4_area_pct",
        "median_stable4_conditional_change_2021_nw",
        "median_stable4_conditional_change_2022_nw",
    ]
    row = {"n_nuts3": int(len(d)), "n_countries": int(d["CNTR_CODE"].nunique())}
    for col in cols:
        if col in country.columns:
            row[f"country_balanced_{col}"] = float(country[col].mean())

    row.update(
        {
            "regional_median_total_2021_change_nw": float(
                d["total_2021_change_nw"].median()
            ),
            "regional_median_new_positive_contribution_nw": float(
                d["new_2021_contribution"].median()
            ),
            "regional_median_lost_positive_contribution_nw": float(
                d["lost_2021_contribution"].median()
            ),
            "regional_median_persistent_positive_change_nw": float(
                d["persistent_2021_change_nw"].median()
            ),
            "regional_median_transient_fraction_of_new_area": float(
                d["transient_fraction_of_new_area"].median()
            ),
            "regional_median_transient_fraction_of_new_contribution": float(
                d["transient_fraction_of_new_contribution"].median()
            ),
            "regional_median_stable4_conditional_change_2021_nw": float(
                d["stable4_conditional_change_2021_nw"].median()
            ),
            "regional_median_stable4_conditional_change_2022_nw": float(
                d["stable4_conditional_change_2022_nw"].median()
            ),
        }
    )
    return pd.DataFrame([row])


def make_plots(
    by_density: pd.DataFrame,
    by_brightness: pd.DataFrame,
    output_dir: Path,
    dpi: int,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    for table, group, order, filename, xlabel in [
        (
            by_density,
            "density_band",
            ["<30", "30-100", "100-300", "300-1000", ">=1000"],
            "pixel_components_2021_by_density.pdf",
            "Fixed population-density band (people km$^{-2}$)",
        ),
        (
            by_brightness,
            "brightness_quintile",
            Q_LABELS,
            "pixel_components_2021_by_brightness.pdf",
            "Fixed 2014-2016 brightness quintile",
        ),
    ]:
        g = table.copy()
        g[group] = pd.Categorical(g[group], categories=order, ordered=True)
        g = g.sort_values(group)
        x = np.arange(len(g))
        width = 0.25
        fig, ax = plt.subplots(figsize=(9.2, 5.9))
        ax.bar(
            x - width,
            g["median_new_positive_contribution_nw"],
            width=width,
            label="new-positive pixels",
        )
        ax.bar(
            x,
            g["median_lost_positive_contribution_nw"],
            width=width,
            label="lost-positive pixels",
        )
        ax.bar(
            x + width,
            g["median_persistent_positive_change_nw"],
            width=width,
            label="persistent-positive pixels",
        )
        ax.axhline(0, linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(g[group].astype(str))
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Median 2020→2021 contribution (nW cm$^{-2}$ sr$^{-1}$)")
        ax.set_title("Direct pixel-support decomposition of the 2021 pulse")
        ax.legend()
        ax.grid(alpha=0.2, axis="y")
        fig.tight_layout()
        fig.savefig(output_dir / filename, dpi=dpi)
        plt.close(fig)

    g = by_density.copy()
    g["density_band"] = pd.Categorical(
        g["density_band"],
        categories=["<30", "30-100", "100-300", "300-1000", ">=1000"],
        ordered=True,
    )
    g = g.sort_values("density_band")
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    ax.bar(
        np.arange(len(g)),
        100.0 * g["median_transient_fraction_of_new_area"],
    )
    ax.set_xticks(np.arange(len(g)))
    ax.set_xticklabels(g["density_band"].astype(str))
    ax.set_xlabel("Fixed population-density band (people km$^{-2}$)")
    ax.set_ylabel("Median fraction of 2021 new-positive area gone in 2022 (%)")
    ax.set_title("Persistence of pixels newly admitted to the 2021 mask")
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(output_dir / "new_2021_pixel_persistence_by_density.pdf", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    x = np.arange(len(g))
    width = 0.36
    ax.bar(
        x - width / 2,
        g["median_stable4_conditional_change_2021_nw"],
        width=width,
        label="2020→2021",
    )
    ax.bar(
        x + width / 2,
        g["median_stable4_conditional_change_2022_nw"],
        width=width,
        label="2021→2022",
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(g["density_band"].astype(str))
    ax.set_xlabel("Fixed population-density band (people km$^{-2}$)")
    ax.set_ylabel(
        "Median radiance change on pixels positive in all 2020-2023 years "
        "(nW cm$^{-2}$ sr$^{-1}$)"
    )
    ax.set_title("Intensity change on temporally stable lit support")
    ax.legend()
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(output_dir / "stable_pixels_2021_2022_by_density.pdf", dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/pipeline.yaml")
    args = ap.parse_args()

    integration_self_test()

    config = load_config(args.config)
    cfg = config["analysis"]["fixed_nuts2016_2021_diagnostic"]

    panel = attach_fixed_brightness(load_locked_panel(config))
    components, manifest = extract_or_load(panel, config)
    regional = build_regional(panel, components)

    expected_n = int(cfg["expected_n_nuts3"])
    if len(regional) != expected_n:
        raise RuntimeError(
            f"Expected {expected_n} regions; found {len(regional)}."
        )

    by_density = summarize_group(regional, "density_band")
    by_brightness = summarize_group(regional, "brightness_quintile")
    country = country_summary(regional)
    overall = overall_summary(regional, country)

    processed = configured_path(config, "processed")
    figures = (
        configured_path(config, "figures")
        / "population_relationship"
        / "fixed_nuts2016_2021_pixel_support"
    )
    processed.mkdir(parents=True, exist_ok=True)

    outputs = {
        "manifest": processed
        / "fixed_nuts2016_2021_pixel_support_manifest.csv",
        "regional": processed
        / "fixed_nuts2016_2021_pixel_support_regional.csv",
        "overall": processed
        / "fixed_nuts2016_2021_pixel_support_overall.csv",
        "country": processed
        / "fixed_nuts2016_2021_pixel_support_by_country.csv",
        "density": processed
        / "fixed_nuts2016_2021_pixel_support_by_density.csv",
        "brightness": processed
        / "fixed_nuts2016_2021_pixel_support_by_brightness.csv",
        "summary": processed
        / "fixed_nuts2016_2021_pixel_support_summary.txt",
    }

    manifest.to_csv(outputs["manifest"], index=False)
    regional.to_csv(outputs["regional"], index=False)
    overall.to_csv(outputs["overall"], index=False)
    country.to_csv(outputs["country"], index=False)
    by_density.to_csv(outputs["density"], index=False)
    by_brightness.to_csv(outputs["brightness"], index=False)

    make_plots(
        by_density,
        by_brightness,
        figures,
        int(cfg["figure_dpi"]),
    )

    o = overall.iloc[0]
    with outputs["summary"].open("w", encoding="utf-8") as h:
        h.write("FIXED NUTS2016 2021 PIXEL-SUPPORT PERSISTENCE AUDIT\n")
        h.write("=" * 80 + "\n\n")
        h.write(
            f"Locked panel: {len(regional)} NUTS3 regions; "
            f"{regional['CNTR_CODE'].nunique()} countries.\n"
        )
        h.write(
            "Direct 2020->2021 pixel-state decomposition of raw "
            "average_masked radiance. Newly positive pixels are followed "
            "into 2022; a stable-support intensity test uses pixels positive "
            "in every year 2020-2023.\n\n"
        )

        h.write("OVERALL SUMMARY\n")
        h.write(overall.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
        h.write("\n\nBY FIXED POPULATION-DENSITY BAND\n")
        h.write(
            by_density.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )
        h.write("\n\nBY FIXED BASELINE BRIGHTNESS\n")
        h.write(
            by_brightness.to_string(index=False, float_format=lambda x: f"{x:.6g}")
        )

        h.write("\n\nDECISION RULE\n")
        h.write(
            "- If a large fraction of 2021 newly positive area disappears in "
            "2022, the mask expansion is predominantly transient rather than "
            "persistent spatial growth of lit area.\n"
        )
        h.write(
            "- If 2021 radiance also rises on pixels positive in every year "
            "2020-2023 and reverses there in 2022, an intensity/product "
            "disturbance exists independently of changing mask support.\n"
        )
        h.write(
            "- If stable-support intensity is near zero while transient new "
            "pixels explain most of the pulse, masking is the dominant "
            "mechanism.\n"
        )
        h.write(
            "- These diagnostics establish mechanism; they do not by "
            "themselves prescribe a correction.\n"
        )

    print("\nFIXED NUTS2016 2021 PIXEL-SUPPORT PERSISTENCE AUDIT")
    print("\nOVERALL SUMMARY")
    print(overall.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nBY FIXED POPULATION-DENSITY BAND")
    print(by_density.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\nWROTE")
    for p in [*outputs.values(), figures]:
        print(p)


if __name__ == "__main__":
    main()
