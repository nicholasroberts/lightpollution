#!/usr/bin/env python3
"""Merge annual NUTS3 VIIRS radiance and population-density tables."""

from __future__ import annotations

import argparse

from project_config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)

    raise NotImplementedError(
        "Merge stage will be enabled after the radiance and population schemas "
        f"are fixed. Target output: {config['outputs']['merged']}."
    )


if __name__ == "__main__":
    main()
