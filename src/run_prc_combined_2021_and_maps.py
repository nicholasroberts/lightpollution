#!/usr/bin/env python3
"""Combined analysis: 2021 corrected PRC and endpoint radiance/PRC graphics.

Run self-tests, then the 2021 spatial overlay using the existing Europe patch,
and produce the two 2014->2024 map pairs + country-faceted PRC/radiance
scatter + joined NUTS3 table. All outputs in separate directories.
No full 2014-2024 raster re-extraction, original satellite data preserved.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
from project_config import PROJECT_ROOT


def step(script,*options):
    call=[sys.executable,str(PROJECT_ROOT/"src"/script),*options]
    print("\nRunning:", " ".join(map(str,call)),flush=True)
    subprocess.run(call,cwd=PROJECT_ROOT,check=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true",
                    help="Validate both stages without input rasters")
    args=ap.parse_args()
    overlay="prc_2021_locked17_overlay_analysis.py"
    figures="plot_prc_radiance_2014_2024.py"
    step(overlay,"--self-test")
    step(figures,"--self-test")
    if args.self_test:
        print("Combined PRC and radiance comparison self-tests: PASS")
        return
    step(overlay)
    step(figures)
    print("\nComplete: candidate 2021 PRC scenario plus four endpoint outputs.")
    print("2021 results: data/processed/prc_2021_locked17_combined/")
    print("2014–2024 maps and country plots: figures/prc_radiance_comparison_locked17/")


if __name__=="__main__":
    main()
