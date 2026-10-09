#!/usr/bin/env python3
"""Validate country-specific VIIRS trends against ISS spectral whitening.

This is deliberately an independent validation rather than an annual LED
penetration covariate. Sánchez de Miguel et al. (2022) published calibrated
ISS RGB mosaics for Europe for 2012-2013 and 2014-2020. Their B/G and G/R
ratios capture the transition from sodium-dominated spectra toward whiter,
bluer lighting.

Workflow:
  1. download the two calibrated ~421 MB RGBA mosaics from Zenodo;
  2. calculate country-level pre/post B/G and G/R ratios from the mosaics;
  3. fit a 2013-2020 NUTS3 mixed model controlling population density,
     GDP per capita PPS, country baseline and NUTS3 random intercept;
  4. estimate each country's residual VIIRS temporal slope;
  5. test whether spectral whitening predicts those residual slopes.

The hypothesis is directional: countries with a stronger shift toward the blue
should, all else equal, show a more negative VIIRS-band temporal trend because
VIIRS DNB is relatively insensitive to the blue peak of white LEDs.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import urllib.request
import warnings
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from scipy.stats import pearsonr, spearmanr
import statsmodels.api as sm
import statsmodels.formula.api as smf

from project_config import configured_path, load_config


ZENODO_ROOT = "https://zenodo.org/records/7677478/files/"
PRE_FILE = "Corr_EUg0RGBA_pre2013_48.tiff"
POST_FILE = "Corr_EUg0RGBA_post2013_46.tif"


def download(url: str, destination: Path, force: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        print(f"Already present: {destination}")
        return

    if shutil.which("wget"):
        command = [
            "wget",
            "--continue",
            "--output-document",
            str(destination),
            url,
        ]
        print("+", " ".join(command), flush=True)
        subprocess.run(command, check=True)
    else:
        print("wget not found; falling back to Python urllib.", flush=True)
        with urllib.request.urlopen(url) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)


def load_country_geometries(config: dict) -> gpd.GeoDataFrame:
    """Use the latest configured NUTS release in which each country occurs."""
    id_field = config["nuts"]["id_field"]
    root = configured_path(config, "nuts_raw")
    template = config["nuts"]["filename_template"]

    latest = {}
    for release in sorted(config["nuts"]["releases"]):
        path = root / template.format(release=release)
        if not path.exists():
            raise FileNotFoundError(
                f"Missing NUTS geometry: {path}. "
                "Run: python src/run_pipeline.py nuts"
            )
        g = gpd.read_file(path)
        if "CNTR_CODE" not in g.columns:
            g["CNTR_CODE"] = g[id_field].astype(str).str[:2]

        for country, sub in g.groupby("CNTR_CODE"):
            geom = sub.geometry.union_all()
            latest[str(country)] = {
                "country": str(country),
                "nuts_release": int(release),
                "geometry": geom,
            }

    return gpd.GeoDataFrame(
        list(latest.values()),
        geometry="geometry",
        crs="EPSG:4326",
    )


def ratio_stats_for_country(
    src: rasterio.io.DatasetReader,
    geometry,
) -> dict:
    geom = gpd.GeoSeries([geometry], crs="EPSG:4326").to_crs(src.crs).iloc[0]

    indexes = [1, 2, 3]
    if src.count >= 4:
        indexes.append(4)

    data, _ = mask(
        src,
        [geom.__geo_interface__],
        crop=True,
        indexes=indexes,
        filled=False,
    )

    r = np.ma.asarray(data[0], dtype=float)
    g = np.ma.asarray(data[1], dtype=float)
    b = np.ma.asarray(data[2], dtype=float)

    base_mask = (
        np.ma.getmaskarray(r)
        | np.ma.getmaskarray(g)
        | np.ma.getmaskarray(b)
        | ~np.isfinite(np.ma.filled(r, np.nan))
        | ~np.isfinite(np.ma.filled(g, np.nan))
        | ~np.isfinite(np.ma.filled(b, np.nan))
        | (np.ma.filled(r, 0) <= 0)
        | (np.ma.filled(g, 0) <= 0)
        | (np.ma.filled(b, 0) < 0)
    )

    if src.count >= 4:
        alpha = np.ma.asarray(data[3], dtype=float)
        base_mask |= (
            np.ma.getmaskarray(alpha)
            | (np.ma.filled(alpha, 0) <= 0)
        )

    rv = np.ma.filled(r, np.nan)
    gv = np.ma.filled(g, np.nan)
    bv = np.ma.filled(b, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        bg = bv / gv
        gr = gv / rv
        rg = rv / gv

    # Apply the published colour-ratio quality bounds used in the paper.
    valid = (
        ~base_mask
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


def extract_country_spectral_change(
    config: dict,
    pre_path: Path,
    post_path: Path,
) -> pd.DataFrame:
    countries = load_country_geometries(config)

    rows = []
    with rasterio.open(pre_path) as pre, rasterio.open(post_path) as post:
        print(
            f"Pre mosaic: {pre.width} x {pre.height}; "
            f"bands={pre.count}; CRS={pre.crs}"
        )
        print(
            f"Post mosaic: {post.width} x {post.height}; "
            f"bands={post.count}; CRS={post.crs}"
        )

        for row in countries.itertuples(index=False):
            try:
                before = ratio_stats_for_country(pre, row.geometry)
                after = ratio_stats_for_country(post, row.geometry)
            except ValueError:
                # No raster overlap.
                continue

            if (
                before["n_valid_pixels"] == 0
                and after["n_valid_pixels"] == 0
            ):
                continue

            record = {
                "CNTR_CODE": row.country,
                "nuts_release_geometry": row.nuts_release,
            }
            for key, value in before.items():
                record[f"pre_{key}"] = value
            for key, value in after.items():
                record[f"post_{key}"] = value

            record["delta_bg_median"] = (
                record["post_bg_median"] - record["pre_bg_median"]
            )
            record["delta_gr_median"] = (
                record["post_gr_median"] - record["pre_gr_median"]
            )
            if record["pre_bg_median"] > 0:
                record["pct_change_bg_median"] = (
                    100
                    * record["delta_bg_median"]
                    / record["pre_bg_median"]
                )
            else:
                record["pct_change_bg_median"] = np.nan

            if record["pre_gr_median"] > 0:
                record["pct_change_gr_median"] = (
                    100
                    * record["delta_gr_median"]
                    / record["pre_gr_median"]
                )
            else:
                record["pct_change_gr_median"] = np.nan

            rows.append(record)
            print(
                f"{row.country}: "
                f"B/G {record['pre_bg_median']:.4f} -> "
                f"{record['post_bg_median']:.4f}; "
                f"G/R {record['pre_gr_median']:.4f} -> "
                f"{record['post_gr_median']:.4f}"
            )

    return pd.DataFrame(rows).sort_values("CNTR_CODE").reset_index(drop=True)


def fit_residual_country_trends(
    config: dict,
    spectral: pd.DataFrame,
) -> tuple[pd.DataFrame, object]:
    """Estimate country temporal trends after population/GDP are controlled."""
    processed = configured_path(config, "processed")
    data = pd.read_parquet(
        processed / config["outputs"]["merged_socioeconomic"]
    ).copy()

    cfg = config["analysis"]["iss_spectral_validation"]
    first = int(cfg["first_year"])
    last = int(cfg["last_year"])
    gdp_metric = config["socioeconomic"]["model_gdp_metric"]
    log_gdp_col = f"log10_{gdp_metric}"

    data = data[
        data["year"].between(first, last, inclusive="both")
        & data["CNTR_CODE"].isin(spectral["CNTR_CODE"])
    ].copy()

    required = [
        "log10_radiance_mean_corrected",
        "log10_population_density",
        log_gdp_col,
    ]
    finite = np.ones(len(data), dtype=bool)
    for col in required:
        finite &= np.isfinite(pd.to_numeric(data[col], errors="coerce"))
    data = data.loc[finite].copy()

    counts = data.groupby("CNTR_CODE")["NUTS_ID"].nunique()
    keep = counts[
        counts >= int(cfg["min_country_regions"])
    ].index
    data = data[data["CNTR_CODE"].isin(keep)].copy()

    repeated = data.groupby("NUTS_ID").size()
    data = data[
        data["NUTS_ID"].isin(repeated[repeated >= 2].index)
    ].copy()

    data["country"] = data["CNTR_CODE"].astype("category")
    data["log_rad"] = data["log10_radiance_mean_corrected"]
    data["log_pop"] = data["log10_population_density"]
    data["log_gdp_pc"] = data[log_gdp_col]

    # Center demand variables. The validation model deliberately contains no
    # country-specific temporal term: the remaining country residual slope is
    # what the ISS spectral transition is being asked to explain.
    data["log_pop"] -= data["log_pop"].mean()
    data["log_gdp_pc"] -= data["log_gdp_pc"].mean()

    formula = (
        "log_rad ~ log_pop + log_gdp_pc + "
        "log_pop:log_gdp_pc + C(country)"
    )
    model = smf.mixedlm(
        formula,
        data=data,
        groups=data["NUTS_ID"],
        re_formula="1",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = model.fit(
            reml=False,
            method="lbfgs",
            maxiter=2000,
            disp=False,
        )

    # Residuals include the fitted NUTS3 random intercept. A random intercept
    # changes level but not temporal slope, which is exactly what we want.
    data["residual"] = result.resid
    midpoint = (first + last) / 2
    data["year_centered_validation"] = data["year"] - midpoint

    rows = []
    for country, sub in data.groupby("CNTR_CODE"):
        if sub["year"].nunique() < 3:
            continue
        X = sm.add_constant(sub["year_centered_validation"])
        fit = sm.OLS(sub["residual"], X).fit()
        rows.append(
            {
                "CNTR_CODE": country,
                "n_nuts3_year_observations": int(len(sub)),
                "n_nuts3_regions": int(sub["NUTS_ID"].nunique()),
                "n_years": int(sub["year"].nunique()),
                "residual_viirs_slope_per_year": float(
                    fit.params["year_centered_validation"]
                ),
                "residual_viirs_slope_se": float(
                    fit.bse["year_centered_validation"]
                ),
                "residual_viirs_slope_p": float(
                    fit.pvalues["year_centered_validation"]
                ),
            }
        )

    return pd.DataFrame(rows), result


def association_table(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    y = data["residual_viirs_slope_per_year"]

    for xcol in ["delta_bg_median", "delta_gr_median"]:
        sub = data[[xcol, "residual_viirs_slope_per_year"]].dropna()
        if len(sub) < 4:
            continue

        x = sub[xcol]
        yy = sub["residual_viirs_slope_per_year"]
        pr, pp = pearsonr(x, yy)
        sr, sp = spearmanr(x, yy)
        X = sm.add_constant(x)
        fit = sm.OLS(yy, X).fit()

        rows.append(
            {
                "spectral_predictor": xcol,
                "n_countries": int(len(sub)),
                "pearson_r": float(pr),
                "pearson_p": float(pp),
                "spearman_rho": float(sr),
                "spearman_p": float(sp),
                "ols_intercept": float(fit.params["const"]),
                "ols_slope": float(fit.params[xcol]),
                "ols_slope_se": float(fit.bse[xcol]),
                "ols_slope_p": float(fit.pvalues[xcol]),
                "ols_r2": float(fit.rsquared),
            }
        )

    complete = data[
        [
            "delta_bg_median",
            "delta_gr_median",
            "residual_viirs_slope_per_year",
        ]
    ].dropna()
    if len(complete) >= 6:
        X = sm.add_constant(complete[["delta_bg_median", "delta_gr_median"]])
        fit = sm.OLS(
            complete["residual_viirs_slope_per_year"],
            X,
        ).fit()
        rows.append(
            {
                "spectral_predictor": "delta_bg_median + delta_gr_median",
                "n_countries": int(len(complete)),
                "pearson_r": np.nan,
                "pearson_p": np.nan,
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
                "ols_intercept": float(fit.params["const"]),
                "ols_slope": np.nan,
                "ols_slope_se": np.nan,
                "ols_slope_p": float(fit.f_pvalue),
                "ols_r2": float(fit.rsquared),
            }
        )

    return pd.DataFrame(rows)


def make_plot(data: pd.DataFrame, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    for xcol, xlabel, filename in [
        (
            "delta_bg_median",
            "Change in ISS B/G ratio (post − pre)",
            "iss_delta_bg_vs_residual_viirs_trend.png",
        ),
        (
            "delta_gr_median",
            "Change in ISS G/R ratio (post − pre)",
            "iss_delta_gr_vs_residual_viirs_trend.png",
        ),
    ]:
        sub = data[[xcol, "residual_viirs_slope_per_year", "CNTR_CODE"]].dropna()
        if len(sub) < 4:
            continue

        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        ax.scatter(
            sub[xcol],
            sub["residual_viirs_slope_per_year"],
            s=28,
        )

        X = sm.add_constant(sub[xcol])
        fit = sm.OLS(sub["residual_viirs_slope_per_year"], X).fit()
        xx = np.linspace(sub[xcol].min(), sub[xcol].max(), 100)
        yy = fit.params["const"] + fit.params[xcol] * xx
        ax.plot(xx, yy)

        for row in sub.itertuples(index=False):
            ax.annotate(
                row.CNTR_CODE,
                (getattr(row, xcol), row.residual_viirs_slope_per_year),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )

        ax.axhline(0, linewidth=0.8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(
            "Residual corrected-VIIRS log10 radiance trend (yr⁻¹)\n"
            "after population density, GDP and country baseline"
        )
        ax.set_title(
            f"ISS spectral shift vs VIIRS residual trend\n"
            f"2013–2020; R² = {fit.rsquared:.3f}, p = {fit.f_pvalue:.3g}"
        )
        fig.tight_layout()
        path = outdir / filename
        fig.savefig(path, dpi=300)
        plt.close(fig)
        print(f"Wrote: {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["iss_spectral_validation"]

    raw_dir = Path(cfg["raw_dir"])
    if not raw_dir.is_absolute():
        raw_dir = Path.cwd() / raw_dir
    raw_dir.mkdir(parents=True, exist_ok=True)

    pre_path = raw_dir / PRE_FILE
    post_path = raw_dir / POST_FILE

    download(
        ZENODO_ROOT + PRE_FILE + "?download=1",
        pre_path,
        force=args.force_download,
    )
    download(
        ZENODO_ROOT + POST_FILE + "?download=1",
        post_path,
        force=args.force_download,
    )

    spectral = extract_country_spectral_change(
        config,
        pre_path,
        post_path,
    )

    processed = configured_path(config, "processed")
    spectral_path = processed / "iss_country_spectral_shift.csv"
    spectral.to_csv(spectral_path, index=False)

    trends, control_model = fit_residual_country_trends(
        config,
        spectral,
    )
    combined = trends.merge(
        spectral,
        on="CNTR_CODE",
        how="inner",
        validate="one_to_one",
    )
    combined_path = processed / "iss_spectral_viirs_country_validation.csv"
    combined.to_csv(combined_path, index=False)

    associations = association_table(combined)
    association_path = processed / "iss_spectral_viirs_associations.csv"
    associations.to_csv(association_path, index=False)

    summary_path = processed / "iss_spectral_validation_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("ISS SPECTRAL-SHIFT / VIIRS TEMPORAL VALIDATION\n")
        handle.write("=" * 72 + "\n\n")
        handle.write(
            "ISS source: Sánchez de Miguel et al. 2022, calibrated Europe "
            "RGB mosaics, Zenodo 10.5281/zenodo.7677478\n"
        )
        handle.write(
            f"Validation interval: {cfg['first_year']}-{cfg['last_year']}\n"
        )
        handle.write(
            f"Countries with spectral + model data: {len(combined)}\n\n"
        )
        handle.write("CONTROL MIXED MODEL\n")
        handle.write(control_model.summary().as_text())
        handle.write("\n\nSPECTRAL ASSOCIATIONS\n")
        handle.write(associations.to_string(index=False))
        handle.write("\n")

    figure_dir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "spectral_validation"
    )
    make_plot(combined, figure_dir)

    print("\nCOUNTRY SPECTRAL SHIFT")
    print(
        spectral[
            [
                "CNTR_CODE",
                "pre_bg_median",
                "post_bg_median",
                "delta_bg_median",
                "pre_gr_median",
                "post_gr_median",
                "delta_gr_median",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.5g}")
    )

    print("\nASSOCIATION WITH RESIDUAL VIIRS TREND")
    print(
        associations.to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )
    print()
    print(f"Wrote: {spectral_path}")
    print(f"Wrote: {combined_path}")
    print(f"Wrote: {association_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
