#!/usr/bin/env python3
"""Merge NUTS3 socioeconomic metrics onto the production analysis table."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from project_config import configured_path, load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")

    base = pd.read_parquet(
        processed / config["outputs"]["merged"]
    ).copy()
    socioeconomic = pd.read_parquet(
        processed / config["outputs"]["socioeconomic_by_nuts3"]
    ).copy()

    id_field = config["nuts"]["id_field"]
    keys = ["year", id_field, "nuts_release"]

    # CNTR_CODE is already present in the production table. Avoid duplicate
    # geography columns while retaining every socioeconomic metric/provenance.
    drop_cols = [
        col for col in ["CNTR_CODE"]
        if col in socioeconomic.columns
    ]
    socioeconomic = socioeconomic.drop(columns=drop_cols)

    merged = base.merge(
        socioeconomic,
        on=keys,
        how="left",
        validate="one_to_one",
        suffixes=("", "_socioeconomic"),
    )

    metric = config["socioeconomic"]["model_gdp_metric"]
    if metric not in merged.columns:
        raise KeyError(f"Configured GDP metric not found after merge: {metric}")

    log_metric = f"log10_{metric}"
    if log_metric not in merged.columns:
        merged[log_metric] = np.nan
        ok = np.isfinite(merged[metric]) & (merged[metric] > 0)
        merged.loc[ok, log_metric] = np.log10(merged.loc[ok, metric])

    output = (
        processed
        / config["outputs"]["merged_socioeconomic"]
    )
    merged.to_parquet(output, index=False)
    merged.to_csv(output.with_suffix(".csv"), index=False)

    summary = (
        merged.groupby("year")
        .agg(
            analysis_rows=(id_field, "size"),
            gdp_rows=(metric, "count"),
            complete_for_max_model=(
                log_metric,
                lambda s: int(
                    (
                        s.notna()
                        & merged.loc[s.index, "log10_population_density"].notna()
                        & merged.loc[
                            s.index,
                            "log10_radiance_mean_corrected",
                        ].notna()
                    ).sum()
                ),
            ),
        )
    )

    print("Socioeconomic merge coverage:")
    print(summary.to_string())
    print()
    print(f"Wrote: {output}")
    print(f"Wrote: {output.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
