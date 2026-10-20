#!/usr/bin/env python3
"""One command: full moth history, approved snow correction for BOTH taxa, three-way maps.

Reuses already-generated spider locked17 history and calibrated radiance table,
so only the new MOTH 2014–2024 EOG extraction is repeated. Original outputs
are preserved. The locked 2021 VIIRS snow correction is applied to both taxa as the standard series.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
from project_config import PROJECT_ROOT


def invoke(script,*args):
    command=[sys.executable,str(PROJECT_ROOT/"src"/script),*args]
    print("\nRUN:", " ".join(command),flush=True)
    subprocess.run(command,cwd=PROJECT_ROOT,check=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true",
                    help="Run lightweight diagnostics without EU raster processing")
    args=ap.parse_args()
    invoke("fit_prc_moth_observation_year.py","--self-test")
    invoke("plot_radiance_spider_moth_2014_2024.py","--self-test")
    invoke("prc_2021_locked_moth16_overlay_analysis.py","--self-test")
    invoke("publish_prc_locked_snow.py","--self-test")
    if args.self_test:
        print("Moth calibration, endpoint and snow-overlay tests: PASS",flush=True)
        return
    invoke("historical_prc_pixel_nuts3.py",
           "--transfer-config","config/prc_transfer_moth.yaml",
           "--output-dir","data/processed/historical_prc_2014_2024_locked_moth16",
           "--figures-dir","figures/historical_prc_2014_2024_locked_moth16")
    invoke("plot_radiance_spider_moth_2014_2024.py")
    invoke("prc_2021_locked17_overlay_analysis.py")
    invoke("prc_2021_locked_moth16_overlay_analysis.py")
    invoke("publish_prc_locked_snow.py")
    print("\nMoth full 2014–2024 extraction, R/spider/moth maps and",
          "2021 LOCKED snow-corrected spider and moth series: COMPLETE.",flush=True)


if __name__=="__main__":
    main()
