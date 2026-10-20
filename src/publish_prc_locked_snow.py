#!/usr/bin/env python3
"""Publish approved 2021 snow-corrected PRC time series for both taxa."""
import argparse
from pathlib import Path
import pandas as pd
import yaml
from project_config import PROJECT_ROOT as ROOT
from historical_prc_pixel_nuts3 import build_panel,europe_table,load_nuts,load_config,plot_region_atlases,plot_europe

def main():
    raise RuntimeError('Publisher implementation pending; do not use uncorrected standard outputs')

if __name__=='__main__':
    main()
