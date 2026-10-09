#!/usr/bin/env python3
"""ISS spectral-shift validation and NUTS3 slope/intercept mechanism analysis.

Sánchez de Miguel et al. (2022) published calibrated ISS RGB mosaics for
Europe for 2012-2013 and 2014-2020. This stage uses those mosaics to ask
whether the spatial pattern of spectral whitening can explain part of the
country-specific change in the VIIRS radiance ~ population-density log-log
relationship.

Key design points
-----------------
1. Restrict ISS analysis to the countries retained in the published paper.
2. Reproduce the published colour-ratio quality filters.
3. Apply a VIIRS-DNB > 0.5 nW cm-2 sr-1 mask to suppress unreliable ratios in
   dark pixels, as in the paper.
4. Extract the spectral change at NUTS3 as well as country level.
5. Test whether spectral change has a population-density gradient within
   countries: this is the mechanism that can alter the observed log-log slope.
6. Compare that spectral-density gradient with the observed 2013 -> 2020
   change in country log-log slopes.
7. Compare country-average spectral shift with change in radiance at a common
   reference population density (log10 density = 2, i.e. 100 people km-2).
8. Fit a NUTS3 change model containing spectral shift x population density and
   ask whether country-specific residual slope heterogeneity remains.

Important temporal approximation
--------------------------------
The original paper used VIIRS imagery chosen to match the dates of individual
ISS images as closely as possible. The public mosaics combine many acquisition
dates. Here, annual 2013 and 2020 VIIRS rasters are used as reproducible
pre/post masks and as endpoint radiance relationships. The post ISS mosaic is
a 2014-2020 composite, so the 2020 comparison is an anchor, not an assertion
that every post-mosaic pixel was acquired in 2020.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import urllib.request
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import geometry_mask, geometry_window
from rasterio.vrt import WarpedVRT
from scipy.stats import linregress, pearsonr, spearmanr
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

from project_config import configured_path, load_config


ZENODO_ROOT = "https://zenodo.org/records/7677478/files/"
PRE_FILE = "Corr_EUg0RGBA_pre2013_48.tiff"
POST_FILE = "Corr_EUg0RGBA_post2013_46.tif"
EXPECTED_MD5 = {
    PRE_FILE: "e0c547d38f42fa7345f6bae26715595f",
    POST_FILE: "6f156899d716f4b01db7d4b26e058dc7",
}
PUBLISHED_EUROPE_MEDIANS = {
    "pre": {"bg": 0.35, "gr": 0.50},
    "post": {"bg": 0.36, "gr": 0.55},
}


def md5sum(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_download(
    url: str,
    destination: Path,
    expected_md5: str,
    force: bool = False,
) -> None:
    """Download atomically, resume partial downloads, and verify Zenodo MD5."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_suffix(destination.suffix + ".part")

    if destination.exists() and not force:
        observed = md5sum(destination)
        if observed == expected_md5:
            print(f"Verified existing file: {destination} ({observed})")
            return
        print(
            f"Existing file has wrong MD5 ({observed}); moving to resumable "
            f"partial file {part.name}."
        )
        destination.replace(part)

    if force and destination.exists():
        destination.unlink()
    if force and part.exists():
        part.unlink()

    if shutil.which("wget"):
        command = [
            "wget",
            "--continue",
            "--output-document",
            str(part),
            url,
        ]
        print("+", " ".join(command), flush=True)
        subprocess.run(command, check=True)
    else:
        if part.exists():
            raise RuntimeError(
                "A partial download exists but wget is unavailable; install "
                "wget or remove the .part file and retry."
            )
        print("wget not found; falling back to Python urllib.", flush=True)
        with urllib.request.urlopen(url) as response, part.open("wb") as out:
            shutil.copyfileobj(response, out)

    observed = md5sum(part)
    if observed != expected_md5:
        raise RuntimeError(
            f"MD5 mismatch for {part}: observed {observed}, "
            f"expected {expected_md5}. The .part file was retained."
        )

    part.replace(destination)
    print(f"Verified download: {destination} ({observed})")


def find_viirs_raster(config: dict, year: int) -> Path:
    root = configured_path(config, "viirs_raw")
    matches = sorted(root.glob(f"*{year}*average_masked*.tif"))
    if not matches:
        raise FileNotFoundError(
            f"No VIIRS average_masked raster found for {year} in {root}."
        )

    if len(matches) > 1:
        # Prefer the newest named product if multiple versions are present.
        preferred = [p for p in matches if "v22" in p.name]
        if preferred:
            matches = preferred
        else:
            preferred = [p for p in matches if "v21" in p.name]
            if preferred:
                matches = preferred

    chosen = sorted(matches)[-1]
    print(f"VIIRS {year} mask/endpoint raster: {chosen}")
    return chosen


def load_nuts2016(config: dict, countries: set[str]) -> gpd.GeoDataFrame:
    release = int(config["analysis"]["iss_spectral_validation"]["nuts_release"])
    path = (
        configured_path(config, "nuts_raw")
        / config["nuts"]["filename_template"].format(release=release)
    )
    if not path.exists():
        raise FileNotFoundError(
            f"Missing NUTS geometry {path}. Run: python src/run_pipeline.py nuts"
        )

    nuts = gpd.read_file(path)
    id_field = config["nuts"]["id_field"]
    if "CNTR_CODE" not in nuts.columns:
        nuts["CNTR_CODE"] = nuts[id_field].astype(str).str[:2]

    nuts = nuts[nuts["CNTR_CODE"].isin(countries)].copy()
    nuts = nuts[[id_field, "CNTR_CODE", "geometry"]].rename(
        columns={id_field: "NUTS_ID"}
    )
    return nuts


def ratio_stats_for_geometry(
    iss: rasterio.io.DatasetReader,
    viirs: rasterio.io.DatasetReader,
    geometry,
    threshold: float,
) -> dict:
    """Calculate paper-filtered ISS ratios inside one polygon."""
    if iss.crs is None:
        raise RuntimeError("ISS mosaic has no CRS.")

    geom = gpd.GeoSeries([geometry], crs="EPSG:4326").to_crs(iss.crs).iloc[0]

    try:
        window = geometry_window(iss, [geom.__geo_interface__])
    except Exception:
        return {
            "n_valid_pixels": 0,
            "bg_median": np.nan,
            "bg_mean": np.nan,
            "gr_median": np.nan,
            "gr_mean": np.nan,
        }

    # geometry_window can touch the raster edge; clip explicitly.
    full = rasterio.windows.Window(0, 0, iss.width, iss.height)
    try:
        window = window.intersection(full)
    except Exception:
        return {
            "n_valid_pixels": 0,
            "bg_median": np.nan,
            "bg_mean": np.nan,
            "gr_median": np.nan,
            "gr_mean": np.nan,
        }

    indexes = [1, 2, 3] + ([4] if iss.count >= 4 else [])
    data = iss.read(indexes, window=window, masked=True).astype("float32")
    v = viirs.read(1, window=window, masked=True).astype("float32")

    if data.shape[1:] != v.shape:
        raise RuntimeError(
            f"ISS/VIIRS aligned-window shape mismatch: "
            f"{data.shape[1:]} vs {v.shape}"
        )

    transform = iss.window_transform(window)
    inside = geometry_mask(
        [geom.__geo_interface__],
        out_shape=v.shape,
        transform=transform,
        invert=True,
    )

    r = np.ma.asarray(data[0], dtype=float)
    g = np.ma.asarray(data[1], dtype=float)
    b = np.ma.asarray(data[2], dtype=float)

    rv = np.ma.filled(r, np.nan)
    gv = np.ma.filled(g, np.nan)
    bv = np.ma.filled(b, np.nan)
    vv = np.ma.filled(v, np.nan)

    invalid = (
        np.ma.getmaskarray(r)
        | np.ma.getmaskarray(g)
        | np.ma.getmaskarray(b)
        | np.ma.getmaskarray(v)
        | ~inside
        | ~np.isfinite(rv)
        | ~np.isfinite(gv)
        | ~np.isfinite(bv)
        | ~np.isfinite(vv)
        | (rv <= 0)
        | (gv <= 0)
        | (bv < 0)
        | (vv <= threshold)
    )

    if iss.count >= 4:
        alpha = np.ma.asarray(data[3], dtype=float)
        av = np.ma.filled(alpha, 0)
        invalid |= np.ma.getmaskarray(alpha) | (av <= 0)

    with np.errstate(divide="ignore", invalid="ignore"):
        bg = bv / gv
        gr = gv / rv
        rg = rv / gv

    # Published quality limits: B/G <= 1.2; G/R <= 1.2; R/G <= 6.
    valid = (
        ~invalid
        & np.isfinite(bg)
        & np.isfinite(gr)
        & np.isfinite(rg)
        & (bg >= 0)
        & (gr >= 0)
        & (bg <= 1.2)
        & (gr <= 1.2)
        & (rg <= 6.0)
    )

    bgv = bg[valid]
    grv = gr[valid]
    if bgv.size == 0:
        return {
            "n_valid_pixels": 0,
            "bg_median": np.nan,
            "bg_mean": np.nan,
            "gr_median": np.nan,
            "gr_mean": np.nan,
        }

    return {
        "n_valid_pixels": int(bgv.size),
        "bg_median": float(np.median(bgv)),
        "bg_mean": float(np.mean(bgv)),
        "gr_median": float(np.median(grv)),
        "gr_mean": float(np.mean(grv)),
    }


def add_change_columns(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    for metric in ("bg_median", "bg_mean", "gr_median", "gr_mean"):
        pre = f"pre_{metric}"
        post = f"post_{metric}"
        if pre in table and post in table:
            table[f"delta_{metric}"] = table[post] - table[pre]
    return table


def extract_country_and_nuts3(
    config: dict,
    pre_path: Path,
    post_path: Path,
    pre_viirs_path: Path,
    post_viirs_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cfg = config["analysis"]["iss_spectral_validation"]
    paper_countries = set(cfg["paper_countries"])
    threshold = float(cfg["viirs_mask_threshold"])

    nuts = load_nuts2016(config, paper_countries)
    present = set(nuts["CNTR_CODE"].unique())
    missing = sorted(paper_countries - present)
    if missing:
        print(
            "Paper countries absent from Eurostat NUTS geometry and therefore "
            f"not extractable here: {missing}"
        )

    countries = nuts.dissolve(by="CNTR_CODE", as_index=False)[
        ["CNTR_CODE", "geometry"]
    ]

    country_rows = []
    nuts_rows = []

    with (
        rasterio.open(pre_path) as pre,
        rasterio.open(post_path) as post,
        rasterio.open(pre_viirs_path) as pre_viirs_src,
        rasterio.open(post_viirs_path) as post_viirs_src,
    ):
        if not (
            pre.crs == post.crs
            and pre.transform == post.transform
            and pre.width == post.width
            and pre.height == post.height
        ):
            raise RuntimeError("Pre/post ISS mosaics do not share one grid.")

        print(
            f"Pre mosaic: {pre.width} x {pre.height}; bands={pre.count}; "
            f"CRS={pre.crs}"
        )
        print(
            f"Post mosaic: {post.width} x {post.height}; bands={post.count}; "
            f"CRS={post.crs}"
        )

        vrt_kwargs = dict(
            crs=pre.crs,
            transform=pre.transform,
            width=pre.width,
            height=pre.height,
            resampling=Resampling.bilinear,
        )

        with (
            WarpedVRT(pre_viirs_src, **vrt_kwargs) as pre_viirs,
            WarpedVRT(post_viirs_src, **vrt_kwargs) as post_viirs,
        ):
            for row in countries.itertuples(index=False):
                before = ratio_stats_for_geometry(
                    pre, pre_viirs, row.geometry, threshold
                )
                after = ratio_stats_for_geometry(
                    post, post_viirs, row.geometry, threshold
                )
                record = {"CNTR_CODE": row.CNTR_CODE}
                record.update({f"pre_{k}": v for k, v in before.items()})
                record.update({f"post_{k}": v for k, v in after.items()})
                country_rows.append(record)
                print(
                    f"{row.CNTR_CODE}: "
                    f"B/G {record['pre_bg_median']:.4f} -> "
                    f"{record['post_bg_median']:.4f}; "
                    f"G/R {record['pre_gr_median']:.4f} -> "
                    f"{record['post_gr_median']:.4f}; "
                    f"pixels {record['pre_n_valid_pixels']}/"
                    f"{record['post_n_valid_pixels']}"
                )

            total = len(nuts)
            for i, row in enumerate(nuts.itertuples(index=False), start=1):
                before = ratio_stats_for_geometry(
                    pre, pre_viirs, row.geometry, threshold
                )
                after = ratio_stats_for_geometry(
                    post, post_viirs, row.geometry, threshold
                )
                record = {
                    "NUTS_ID": row.NUTS_ID,
                    "CNTR_CODE": row.CNTR_CODE,
                }
                record.update({f"pre_{k}": v for k, v in before.items()})
                record.update({f"post_{k}": v for k, v in after.items()})
                nuts_rows.append(record)
                if i % 100 == 0 or i == total:
                    print(f"NUTS3 spectral extraction: {i}/{total}")

            europe_geom = countries.geometry.union_all()
            sanity_rows = []
            for label, iss, vv in [
                ("pre", pre, pre_viirs),
                ("post", post, post_viirs),
            ]:
                stats = ratio_stats_for_geometry(
                    iss, vv, europe_geom, threshold
                )
                expected = PUBLISHED_EUROPE_MEDIANS[label]
                sanity_rows.append(
                    {
                        "period": label,
                        "n_valid_pixels": stats["n_valid_pixels"],
                        "observed_bg_median": stats["bg_median"],
                        "published_bg_median": expected["bg"],
                        "observed_gr_median": stats["gr_median"],
                        "published_gr_median": expected["gr"],
                    }
                )

    country = add_change_columns(pd.DataFrame(country_rows))
    nuts3 = add_change_columns(pd.DataFrame(nuts_rows))
    sanity = pd.DataFrame(sanity_rows)

    country["source_doi"] = "10.1126/sciadv.abl6891"
    country["source_dataset_doi"] = "10.5281/zenodo.7677478"
    nuts3["source_doi"] = "10.1126/sciadv.abl6891"
    nuts3["source_dataset_doi"] = "10.5281/zenodo.7677478"

    return country, nuts3, sanity


def prepare_endpoint_data(
    config: dict,
    nuts3: pd.DataFrame,
) -> pd.DataFrame:
    cfg = config["analysis"]["iss_spectral_validation"]
    processed = configured_path(config, "processed")
    data = pd.read_parquet(
        processed / config["outputs"]["merged_socioeconomic"]
    ).copy()

    pre_year = int(cfg["endpoint_pre_year"])
    post_year = int(cfg["endpoint_post_year"])
    gdp_metric = config["socioeconomic"]["model_gdp_metric"]
    gdp_col = f"log10_{gdp_metric}"
    rad_col = "log10_radiance_mean_corrected"
    pop_col = "log10_population_density"

    keep_cols = [
        "NUTS_ID",
        "CNTR_CODE",
        "year",
        rad_col,
        pop_col,
        gdp_col,
    ]
    data = data[keep_cols].copy()

    before = data[data["year"] == pre_year].drop(columns="year").rename(
        columns={
            rad_col: "log_rad_pre",
            pop_col: "log_pop_pre",
            gdp_col: "log_gdp_pc_pre",
        }
    )
    after = data[data["year"] == post_year].drop(columns="year").rename(
        columns={
            rad_col: "log_rad_post",
            pop_col: "log_pop_post",
            gdp_col: "log_gdp_pc_post",
        }
    )

    paired = before.merge(
        after,
        on=["NUTS_ID", "CNTR_CODE"],
        how="inner",
        validate="one_to_one",
    )
    paired = paired.merge(
        nuts3,
        on=["NUTS_ID", "CNTR_CODE"],
        how="inner",
        validate="one_to_one",
    )

    min_pixels = int(cfg["min_valid_pixels_nuts3"])
    paired = paired[
        (paired["pre_n_valid_pixels"] >= min_pixels)
        & (paired["post_n_valid_pixels"] >= min_pixels)
    ].copy()

    paired["delta_log_rad"] = paired["log_rad_post"] - paired["log_rad_pre"]
    paired["delta_log_pop"] = paired["log_pop_post"] - paired["log_pop_pre"]
    paired["delta_log_gdp_pc"] = (
        paired["log_gdp_pc_post"] - paired["log_gdp_pc_pre"]
    )

    paired["log_pop_within_pre"] = (
        paired["log_pop_pre"]
        - paired.groupby("CNTR_CODE")["log_pop_pre"].transform("mean")
    )

    for metric in ("delta_bg_median", "delta_gr_median"):
        paired[f"{metric}_within"] = (
            paired[metric]
            - paired.groupby("CNTR_CODE")[metric].transform("mean")
        )

    return paired


def spectral_density_models(
    paired: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    model_rows = []
    gradient_rows = []

    for metric in ("delta_bg_median", "delta_gr_median"):
        data = paired[
            np.isfinite(paired[metric])
            & np.isfinite(paired["log_pop_within_pre"])
        ].copy()

        counts = data.groupby("CNTR_CODE")["NUTS_ID"].nunique()
        keep = counts[counts >= 3].index
        data = data[data["CNTR_CODE"].isin(keep)].copy()

        common_formula = (
            f"{metric} ~ C(CNTR_CODE) + log_pop_within_pre"
        )
        hetero_formula = (
            f"{metric} ~ C(CNTR_CODE) * log_pop_within_pre"
        )
        common_plain = smf.ols(common_formula, data=data).fit()
        hetero_plain = smf.ols(hetero_formula, data=data).fit()
        common_robust = smf.ols(common_formula, data=data).fit(cov_type="HC3")
        comparison = anova_lm(common_plain, hetero_plain)

        model_rows.append(
            {
                "spectral_metric": metric,
                "n_nuts3": len(data),
                "n_countries": data["CNTR_CODE"].nunique(),
                "common_within_country_density_slope": float(
                    common_robust.params["log_pop_within_pre"]
                ),
                "common_slope_se_hc3": float(
                    common_robust.bse["log_pop_within_pre"]
                ),
                "common_slope_p_hc3": float(
                    common_robust.pvalues["log_pop_within_pre"]
                ),
                "common_model_r2": float(common_plain.rsquared),
                "country_specific_gradient_F": float(
                    comparison.iloc[-1]["F"]
                ),
                "country_specific_gradient_p": float(
                    comparison.iloc[-1]["Pr(>F)"]
                ),
                "heterogeneous_model_r2": float(hetero_plain.rsquared),
            }
        )

        for country, sub in data.groupby("CNTR_CODE"):
            if len(sub) < 3 or sub["log_pop_pre"].nunique() < 2:
                continue
            fit = linregress(sub["log_pop_pre"], sub[metric])
            gradient_rows.append(
                {
                    "CNTR_CODE": country,
                    "spectral_metric": metric,
                    "n_nuts3": len(sub),
                    "spectral_change_vs_logpop_slope": float(fit.slope),
                    "spectral_change_vs_logpop_se": float(fit.stderr),
                    "spectral_change_vs_logpop_p": float(fit.pvalue),
                    "spectral_change_vs_logpop_r2": float(fit.rvalue**2),
                }
            )

    return pd.DataFrame(model_rows), pd.DataFrame(gradient_rows)


def country_loglog_changes(
    paired: pd.DataFrame,
    reference_log_pop: float,
) -> pd.DataFrame:
    rows = []
    for country, sub in paired.groupby("CNTR_CODE"):
        good = (
            np.isfinite(sub["log_pop_pre"])
            & np.isfinite(sub["log_pop_post"])
            & np.isfinite(sub["log_rad_pre"])
            & np.isfinite(sub["log_rad_post"])
        )
        sub = sub.loc[good]
        if len(sub) < 3:
            continue
        if sub["log_pop_pre"].nunique() < 2 or sub["log_pop_post"].nunique() < 2:
            continue

        pre = linregress(sub["log_pop_pre"], sub["log_rad_pre"])
        post = linregress(sub["log_pop_post"], sub["log_rad_post"])

        pre_ref = pre.intercept + pre.slope * reference_log_pop
        post_ref = post.intercept + post.slope * reference_log_pop

        rows.append(
            {
                "CNTR_CODE": country,
                "n_nuts3": len(sub),
                "slope_pre": float(pre.slope),
                "slope_post": float(post.slope),
                "delta_slope": float(post.slope - pre.slope),
                "intercept_pre": float(pre.intercept),
                "intercept_post": float(post.intercept),
                "delta_intercept": float(post.intercept - pre.intercept),
                "lograd_at_reference_pre": float(pre_ref),
                "lograd_at_reference_post": float(post_ref),
                "delta_lograd_at_reference": float(post_ref - pre_ref),
                "reference_log10_population_density": reference_log_pop,
            }
        )

    return pd.DataFrame(rows)


def association_row(
    data: pd.DataFrame,
    xcol: str,
    ycol: str,
    label: str,
) -> dict | None:
    sub = data[[xcol, ycol]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(sub) < 4 or sub[xcol].nunique() < 2:
        return None

    pear = pearsonr(sub[xcol], sub[ycol])
    spear = spearmanr(sub[xcol], sub[ycol])
    fit = linregress(sub[xcol], sub[ycol])
    return {
        "association": label,
        "x": xcol,
        "y": ycol,
        "n_countries": len(sub),
        "pearson_r": float(pear.statistic),
        "pearson_p": float(pear.pvalue),
        "spearman_rho": float(spear.statistic),
        "spearman_p": float(spear.pvalue),
        "ols_intercept": float(fit.intercept),
        "ols_slope": float(fit.slope),
        "ols_slope_se": float(fit.stderr),
        "ols_slope_p": float(fit.pvalue),
        "ols_r2": float(fit.rvalue**2),
    }


def slope_intercept_associations(
    country_change: pd.DataFrame,
    country_spectral: pd.DataFrame,
    gradients: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for metric, country_metric in [
        ("delta_bg_median", "delta_bg_median"),
        ("delta_gr_median", "delta_gr_median"),
    ]:
        grad = gradients[
            gradients["spectral_metric"] == metric
        ][
            ["CNTR_CODE", "spectral_change_vs_logpop_slope"]
        ]
        slope_data = country_change.merge(
            grad, on="CNTR_CODE", how="inner"
        )
        row = association_row(
            slope_data,
            "spectral_change_vs_logpop_slope",
            "delta_slope",
            f"{metric}: population-density spectral gradient vs log-log slope change",
        )
        if row:
            rows.append(row)

        level_data = country_change.merge(
            country_spectral[["CNTR_CODE", country_metric]],
            on="CNTR_CODE",
            how="inner",
        )
        row = association_row(
            level_data,
            country_metric,
            "delta_lograd_at_reference",
            f"{metric}: country median spectral shift vs radiance-level change",
        )
        if row:
            rows.append(row)

        row = association_row(
            level_data,
            country_metric,
            "delta_intercept",
            f"{metric}: country median spectral shift vs raw intercept change",
        )
        if row:
            rows.append(row)

    return pd.DataFrame(rows)


def nuts3_change_models(paired: pd.DataFrame) -> pd.DataFrame:
    """Test spectral main effect and spectral x density effect on VIIRS change.

    The final nested comparison asks whether country-specific density slopes
    remain necessary after spectral change is included. If not, the spectral
    transition has accounted for the between-country slope heterogeneity.
    """
    rows = []

    for metric in ("delta_bg_median", "delta_gr_median"):
        spectral_within = f"{metric}_within"

        for model_name, include_gdp in [
            ("population_only", False),
            ("plus_gdp_change", True),
        ]:
            required = [
                "delta_log_rad",
                "delta_log_pop",
                "log_pop_within_pre",
                spectral_within,
            ]
            if include_gdp:
                required.append("delta_log_gdp_pc")

            data = paired.copy()
            finite = np.ones(len(data), dtype=bool)
            for col in required:
                finite &= np.isfinite(pd.to_numeric(data[col], errors="coerce"))
            data = data.loc[finite].copy()

            counts = data.groupby("CNTR_CODE")["NUTS_ID"].nunique()
            keep = counts[counts >= 3].index
            data = data[data["CNTR_CODE"].isin(keep)].copy()

            if len(data) < 20:
                continue

            controls = "delta_log_pop"
            if include_gdp:
                controls += " + delta_log_gdp_pc"

            base_formula = (
                f"delta_log_rad ~ C(CNTR_CODE) + {controls} + "
                f"log_pop_within_pre + {spectral_within} + "
                f"log_pop_within_pre:{spectral_within}"
            )
            hetero_formula = (
                base_formula + " + log_pop_within_pre:C(CNTR_CODE)"
            )

            base_plain = smf.ols(base_formula, data=data).fit()
            base_robust = smf.ols(base_formula, data=data).fit(cov_type="HC3")
            hetero_plain = smf.ols(hetero_formula, data=data).fit()
            comparison = anova_lm(base_plain, hetero_plain)

            interaction = f"log_pop_within_pre:{spectral_within}"
            if interaction not in base_robust.params:
                interaction = f"{spectral_within}:log_pop_within_pre"

            rows.append(
                {
                    "spectral_metric": metric,
                    "model": model_name,
                    "n_nuts3": len(data),
                    "n_countries": data["CNTR_CODE"].nunique(),
                    "spectral_main_effect": float(
                        base_robust.params[spectral_within]
                    ),
                    "spectral_main_se_hc3": float(
                        base_robust.bse[spectral_within]
                    ),
                    "spectral_main_p_hc3": float(
                        base_robust.pvalues[spectral_within]
                    ),
                    "spectral_x_density_effect": float(
                        base_robust.params[interaction]
                    ),
                    "spectral_x_density_se_hc3": float(
                        base_robust.bse[interaction]
                    ),
                    "spectral_x_density_p_hc3": float(
                        base_robust.pvalues[interaction]
                    ),
                    "base_r2": float(base_plain.rsquared),
                    "residual_country_density_F": float(
                        comparison.iloc[-1]["F"]
                    ),
                    "residual_country_density_p": float(
                        comparison.iloc[-1]["Pr(>F)"]
                    ),
                    "heterogeneous_r2": float(hetero_plain.rsquared),
                }
            )

    return pd.DataFrame(rows)


def plot_country_gradient_vs_slope(
    country_change: pd.DataFrame,
    gradients: pd.DataFrame,
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)

    for metric, short in [
        ("delta_bg_median", "bg"),
        ("delta_gr_median", "gr"),
    ]:
        grad = gradients[
            gradients["spectral_metric"] == metric
        ][
            ["CNTR_CODE", "spectral_change_vs_logpop_slope"]
        ]
        data = country_change.merge(grad, on="CNTR_CODE", how="inner").dropna(
            subset=["spectral_change_vs_logpop_slope", "delta_slope"]
        )
        if len(data) < 4:
            continue

        fit = linregress(
            data["spectral_change_vs_logpop_slope"],
            data["delta_slope"],
        )
        xx = np.linspace(
            data["spectral_change_vs_logpop_slope"].min(),
            data["spectral_change_vs_logpop_slope"].max(),
            100,
        )
        yy = fit.intercept + fit.slope * xx

        fig, ax = plt.subplots(figsize=(7.4, 5.4))
        ax.scatter(
            data["spectral_change_vs_logpop_slope"],
            data["delta_slope"],
            s=35,
        )
        ax.plot(xx, yy)
        ax.axhline(0, linewidth=0.8)
        ax.axvline(0, linewidth=0.8)

        for row in data.itertuples(index=False):
            ax.annotate(
                row.CNTR_CODE,
                (
                    row.spectral_change_vs_logpop_slope,
                    row.delta_slope,
                ),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )

        ax.set_xlabel(
            f"Within-country population-density gradient of {metric}"
        )
        ax.set_ylabel(
            "Change in country log-log VIIRS slope, 2013 to 2020"
        )
        ax.set_title(
            f"Spectral rollout gradient vs VIIRS slope change\n"
            f"r = {fit.rvalue:.3f}, p = {fit.pvalue:.3g}"
        )
        fig.tight_layout()
        out = figure_dir / f"iss_{short}_density_gradient_vs_slope_change.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        print(f"Wrote: {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["iss_spectral_validation"]

    raw_dir = Path(cfg["raw_dir"])
    if not raw_dir.is_absolute():
        raw_dir = Path(__file__).resolve().parents[1] / raw_dir
    raw_dir.mkdir(parents=True, exist_ok=True)

    pre_path = raw_dir / PRE_FILE
    post_path = raw_dir / POST_FILE

    ensure_download(
        ZENODO_ROOT + PRE_FILE + "?download=1",
        pre_path,
        EXPECTED_MD5[PRE_FILE],
        force=args.force_download,
    )
    ensure_download(
        ZENODO_ROOT + POST_FILE + "?download=1",
        post_path,
        EXPECTED_MD5[POST_FILE],
        force=args.force_download,
    )

    pre_mask_year = int(cfg["viirs_mask_pre_year"])
    post_mask_year = int(cfg["viirs_mask_post_year"])
    pre_viirs = find_viirs_raster(config, pre_mask_year)
    post_viirs = find_viirs_raster(config, post_mask_year)

    country, nuts3, sanity = extract_country_and_nuts3(
        config,
        pre_path,
        post_path,
        pre_viirs,
        post_viirs,
    )

    processed = configured_path(config, "processed")
    country_path = processed / "iss_country_spectral_shift.csv"
    nuts3_path = processed / "iss_nuts3_spectral_shift.csv"
    sanity_path = processed / "iss_spectral_mask_validation.csv"

    country.to_csv(country_path, index=False)
    nuts3.to_csv(nuts3_path, index=False)
    sanity.to_csv(sanity_path, index=False)

    paired = prepare_endpoint_data(config, nuts3)
    paired_path = processed / "iss_nuts3_spectral_viirs_endpoints.csv"
    paired.to_csv(paired_path, index=False)

    density_models, gradients = spectral_density_models(paired)
    density_models_path = processed / "iss_spectral_density_models.csv"
    gradients_path = processed / "iss_country_spectral_density_gradients.csv"
    density_models.to_csv(density_models_path, index=False)
    gradients.to_csv(gradients_path, index=False)

    reference_log_pop = float(cfg["reference_log10_population_density"])
    country_change = country_loglog_changes(paired, reference_log_pop)
    country_change_path = processed / "iss_country_loglog_change_2013_2020.csv"
    country_change.to_csv(country_change_path, index=False)

    associations = slope_intercept_associations(
        country_change,
        country,
        gradients,
    )
    associations_path = (
        processed / "iss_spectral_vs_loglog_change_associations.csv"
    )
    associations.to_csv(associations_path, index=False)

    change_models = nuts3_change_models(paired)
    change_models_path = processed / "iss_nuts3_viirs_change_models.csv"
    change_models.to_csv(change_models_path, index=False)

    figure_dir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "spectral_validation"
    )
    plot_country_gradient_vs_slope(
        country_change,
        gradients,
        figure_dir,
    )

    summary_path = processed / "iss_spectral_validation_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("ISS SPECTRAL-SHIFT / VIIRS SLOPE-INTERCEPT VALIDATION\n")
        handle.write("=" * 78 + "\n\n")
        handle.write(
            "ISS source: Sánchez de Miguel et al. 2022, Science Advances "
            "8:eabl6891; Zenodo 10.5281/zenodo.7677478\n"
        )
        handle.write(
            f"Paper VIIRS mask threshold: > {cfg['viirs_mask_threshold']} "
            "nW cm-2 sr-1\n"
        )
        handle.write(
            f"Reproducible mask approximation: annual VIIRS "
            f"{pre_mask_year} for pre mosaic and {post_mask_year} for post "
            "mosaic\n"
        )
        handle.write(
            f"NUTS release: {cfg['nuts_release']}; "
            f"minimum valid ISS pixels per NUTS3/period: "
            f"{cfg['min_valid_pixels_nuts3']}\n"
        )
        handle.write(
            f"Endpoint relationship: {cfg['endpoint_pre_year']} -> "
            f"{cfg['endpoint_post_year']}\n\n"
        )

        handle.write("EUROPE-WIDE MASK SANITY CHECK\n")
        handle.write(sanity.to_string(index=False))
        handle.write("\n\nSPECTRAL CHANGE ~ WITHIN-COUNTRY POPULATION DENSITY\n")
        handle.write(density_models.to_string(index=False))
        handle.write("\n\nCOUNTRY SPECTRAL-DENSITY GRADIENTS VS SLOPE/LEVEL CHANGE\n")
        handle.write(associations.to_string(index=False))
        handle.write("\n\nNUTS3 VIIRS CHANGE MODELS\n")
        handle.write(change_models.to_string(index=False))
        handle.write("\n\nINTERPRETATION GUIDE\n")
        handle.write(
            "- A spectral-change x population-density effect tests whether "
            "spatially structured spectral conversion can change the observed "
            "VIIRS log-log slope.\n"
        )
        handle.write(
            "- The residual_country_density_p test asks whether country-"
            "specific density slopes are still needed after spectral change "
            "is included. A large p is the desired simplification result.\n"
        )
        handle.write(
            "- Country median spectral shift is compared with change in "
            "radiance at log10(population density)=2 (100 people km-2) as a "
            "more interpretable level/intercept diagnostic.\n"
        )
        handle.write(
            "- The post ISS mosaic is a 2014-2020 composite; using 2020 as "
            "the endpoint is a transparent anchor rather than exact temporal "
            "matching for every pixel.\n"
        )

    print("\nEUROPE-WIDE MASK VALIDATION")
    print(
        sanity.to_string(
            index=False,
            float_format=lambda x: f"{x:.5g}",
        )
    )
    print("\nSPECTRAL CHANGE ~ WITHIN-COUNTRY POPULATION DENSITY")
    print(
        density_models.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nSPECTRAL GRADIENT / SLOPE-INTERCEPT ASSOCIATIONS")
    print(
        associations.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print("\nNUTS3 VIIRS CHANGE MODELS")
    print(
        change_models.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nWROTE")
    for path in [
        country_path,
        nuts3_path,
        sanity_path,
        paired_path,
        density_models_path,
        gradients_path,
        country_change_path,
        associations_path,
        change_models_path,
        summary_path,
    ]:
        print(path)


if __name__ == "__main__":
    main()
