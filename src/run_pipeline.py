#!/usr/bin/env python3
"""Top-level entry point for the European VIIRS/NUTS3 analysis pipeline."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from project_config import PROJECT_ROOT


STAGES = {
    "nuts": PROJECT_ROOT / "src" / "fetch_nuts_boundaries.py",
    "validate": PROJECT_ROOT / "src" / "validate_inputs.py",
    "weighting-test": PROJECT_ROOT / "src" / "validate_zonal_weighting.py",
    "radiance": PROJECT_ROOT / "src" / "extract_nuts3_radiance.py",
    "population": PROJECT_ROOT / "src" / "fetch_eurostat_population.py",
    "population-fallback-test": PROJECT_ROOT / "src" / "validate_population_fallback.py",
    "population-recover": PROJECT_ROOT / "src" / "recover_population_density.py",
    "merge": PROJECT_ROOT / "src" / "merge_radiance_population.py",
    "fit": PROJECT_ROOT / "src" / "fit_population_radiance.py",
    "benchmark": PROJECT_ROOT / "src" / "compare_legacy_fits.py",
    "trends": PROJECT_ROOT / "src" / "plot_fit_trends.py",
    "calibration-offset-test": PROJECT_ROOT / "src" / "test_calibration_offset.py",
    "calibration-country-test": PROJECT_ROOT / "src" / "test_country_calibration_correction.py",
    "calibration-sensitivity": PROJECT_ROOT / "src" / "calibration_sensitivity.py",
}


def run_script(script: Path, config: str) -> None:
    command = [sys.executable, str(script), "--config", config]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=[*STAGES, "all"])
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    if args.stage == "all":
        order = ["nuts", "validate", "radiance", "population", "merge"]
    else:
        order = [args.stage]

    for stage in order:
        run_script(STAGES[stage], args.config)


if __name__ == "__main__":
    main()
