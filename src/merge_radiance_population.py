#!/usr/bin/env python3
"""Merge annual NUTS3 VIIRS radiance and Eurostat population density.

The merge is intentionally an inner join on year + NUTS_ID + nuts_release.
Radiance comes from the production-calibrated table, which retains both raw
and corrected values. Unmatched or missing regions simply do not contribute
a point to that year's country-level relationship.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from project_config import configured_path, load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    processed = configured_path(config, "processed")

    radiance_path = processed / config["outputs"]["radiance_by_nuts3_calibrated"]
    population_path = processed / config["outputs"]["population_by_nuts3"]

    radiance = pd.read_parquet(radiance_path)
    population = pd.read_parquet(population_path)

    id_field = config["nuts"]["id_field"]

    population = population[
        population["matches_nuts_release"]
        & population["population_density"].notna()
    ].copy()

    merged = radiance.merge(
        population,
        on=["year", id_field, "nuts_release"],
        how="inner",
        validate="one_to_one",
        suffixes=("", "_population"),
    )

    # The corrected value is the production radiance metric. Raw radiance is
    # retained alongside it for provenance and diagnostic comparisons.
    merged["radiance_mean"] = merged["radiance_mean_corrected"]

    # Log-ready fields. Non-positive radiance cannot be represented on a
    # logarithmic axis; retain the linear values and mark log10 as missing.
    merged["log10_radiance_mean_raw"] = np.nan
    raw_positive = merged["radiance_mean_raw"] > 0
    merged.loc[raw_positive, "log10_radiance_mean_raw"] = np.log10(
        merged.loc[raw_positive, "radiance_mean_raw"]
    )

    merged["log10_radiance_mean_corrected"] = np.nan
    corrected_positive = merged["radiance_mean_corrected"] > 0
    merged.loc[corrected_positive, "log10_radiance_mean_corrected"] = np.log10(
        merged.loc[corrected_positive, "radiance_mean_corrected"]
    )

    # Backwards-compatible alias used by downstream analysis: this is now the
    # corrected production radiance.
    merged["log10_radiance_mean"] = merged["log10_radiance_mean_corrected"]
    merged["log10_population_density"] = np.nan
    pop_positive = merged["population_density"] > 0
    merged.loc[pop_positive, "log10_population_density"] = np.log10(
        merged.loc[pop_positive, "population_density"]
    )

    merged = merged.sort_values(["year", "CNTR_CODE", id_field]).reset_index(drop=True)

    parquet_path = processed / config["outputs"]["merged"]
    csv_path = processed / config["outputs"]["csv_export"]

    merged.to_parquet(parquet_path, index=False)
    merged.to_csv(csv_path, index=False)

    counts = (
        merged.groupby(["year", "CNTR_CODE"])
        .size()
        .rename("n_pairs")
        .reset_index()
    )

    print(f"Rows written: {len(merged):,}")
    print(f"Years: {sorted(merged['year'].unique().tolist())}")
    print()
    print("Paired NUTS3 observations by country/year:")
    print(counts.to_string(index=False))
    print()
    print(f"Wrote: {parquet_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
