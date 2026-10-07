#!/usr/bin/env python3
"""Extract annual VIIRS radiance statistics for fixed NUTS3 regions.

The module boundary is intentional. Before implementing production extraction
we will validate the exact polygon-overlap and physical pixel-area weighting
method on representative regions and compare it with simpler alternatives.
"""

from __future__ import annotations

import argparse

from project_config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    method = config["analysis"]["area_weighting"]["method"]

    raise NotImplementedError(
        "NUTS3 radiance extraction is intentionally not enabled yet. "
        "Validate the polygon-overlap and pixel-area weighting method first. "
        f"Configured method: {method!r}."
    )


if __name__ == "__main__":
    main()
