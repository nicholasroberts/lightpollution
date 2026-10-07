#!/usr/bin/env python3
"""Fetch and normalise Eurostat NUTS3 population-density data.

Implementation is deliberately separated from radiance extraction. The target
dataset is configured as DEMO_R_D3DENS and the eventual output will contain
one row per year and NUTS3 identifier.
"""

from __future__ import annotations

import argparse

from project_config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    dataset = config["population"]["dataset"]

    raise NotImplementedError(
        f"Eurostat download stage not enabled yet. Configured dataset: {dataset}."
    )


if __name__ == "__main__":
    main()
