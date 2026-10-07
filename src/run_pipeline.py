#!/usr/bin/env python3
"""Top-level entry point for the European VIIRS/NUTS3 analysis pipeline.

At this stage the runner exposes the intended workflow and implements input
validation. Scientific processing stages are added independently so they can
be tested before being enabled here.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from project_config import PROJECT_ROOT


STAGES = {
    "validate": PROJECT_ROOT / "src" / "validate_inputs.py",
    "radiance": PROJECT_ROOT / "src" / "extract_nuts3_radiance.py",
    "population": PROJECT_ROOT / "src" / "fetch_eurostat_population.py",
    "merge": PROJECT_ROOT / "src" / "merge_radiance_population.py",
}


def run_script(script: Path, config: str) -> None:
    command = [sys.executable, str(script), "--config", config]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage",
        choices=[*STAGES, "all"],
        help="Pipeline stage to run.",
    )
    parser.add_argument(
        "--config",
        default="config/pipeline.yaml",
    )
    args = parser.parse_args()

    if args.stage == "all":
        order = ["validate", "radiance", "population", "merge"]
    else:
        order = [args.stage]

    for stage in order:
        run_script(STAGES[stage], args.config)


if __name__ == "__main__":
    main()
