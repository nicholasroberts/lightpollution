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
    "calibration": PROJECT_ROOT / "src" / "apply_radiance_calibration.py",
    "population": PROJECT_ROOT / "src" / "fetch_eurostat_population.py",
    "population-fallback-test": PROJECT_ROOT / "src" / "validate_population_fallback.py",
    "population-recover": PROJECT_ROOT / "src" / "recover_population_density.py",
    "merge": PROJECT_ROOT / "src" / "merge_radiance_population.py",
    "fit": PROJECT_ROOT / "src" / "fit_population_radiance.py",
    "benchmark": PROJECT_ROOT / "src" / "compare_legacy_fits.py",
    "trends": PROJECT_ROOT / "src" / "plot_fit_trends.py",
    "country-facets": PROJECT_ROOT / "src" / "plot_country_facets.py",
    "mixed-model": PROJECT_ROOT / "src" / "mixed_effects_analysis.py",
    "europe-maps": PROJECT_ROOT / "src" / "plot_europe_radiance_maps.py",
    "socioeconomic": PROJECT_ROOT / "src" / "fetch_eurostat_socioeconomic.py",
    "socioeconomic-merge": PROJECT_ROOT / "src" / "merge_socioeconomic.py",
    "socioeconomic-model": PROJECT_ROOT / "src" / "simplify_socioeconomic_model.py",
    "iss-spectral-validation": PROJECT_ROOT / "src" / "validate_iss_spectral_shift.py",
    "iss-led-switch-density-test": PROJECT_ROOT / "src" / "test_iss_led_switch_density.py",
    "iss-led-switch-annual-slopes": PROJECT_ROOT / "src" / "test_iss_led_switch_annual_slopes.py",
    "led-spectral-plateau": PROJECT_ROOT / "src" / "calculate_led_spectral_plateau.py",
    "density-temporal-change": PROJECT_ROOT / "src" / "analyze_density_temporal_change.py",
    "density-country-decomposition": PROJECT_ROOT / "src" / "decompose_country_structure.py",
    "density-demographic-surface": PROJECT_ROOT / "src" / "fit_density_demographic_surface.py",
    "future-radiance-response": PROJECT_ROOT / "src" / "build_future_radiance_response.py",
    "future-radiance-history-validation": PROJECT_ROOT / "src" / "validate_future_radiance_history.py",
    "future-radiance-exclude-2013": PROJECT_ROOT / "src" / "test_exclude_2013_density_trends.py",
    "future-radiance-fixed-2014-2024": PROJECT_ROOT / "src" / "validate_fixed_nuts2016_2014_2024.py",
    "future-radiance-2021-diagnostic": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_2021.py",
    "future-radiance-2021-reversal": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_2021_reversal.py",
    "future-radiance-product-layer-audit": PROJECT_ROOT / "src" / "diagnose_viirs_annual_product_layers.py",
    "future-radiance-absolute-offset-audit": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_absolute_offsets.py",
    "future-radiance-mask-occupancy-audit": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_mask_occupancy.py",
    "future-radiance-mask-intensity-decomposition": PROJECT_ROOT / "src" / "decompose_fixed_nuts2016_mask_intensity.py",
    "future-radiance-2021-pixel-support": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_2021_pixel_support.py",
    "future-radiance-july-vcm-vcmsl-pilot": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_july_vcm_vcmsl.py",
    "future-radiance-download-monthly-2021": PROJECT_ROOT / "src" / "fetch_eog_monthly_vcmslcfg_2021.py",
    "future-radiance-download-monthly-vcmcfg": PROJECT_ROOT / "src" / "fetch_eog_monthly_vcmcfg_2020_2022.py",
    "future-radiance-monthly-2021-localization": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_monthly_2021_localization.py",
    "future-radiance-monthly-product-atlas": PROJECT_ROOT / "src" / "diagnose_fixed_nuts2016_monthly_product_atlas.py",
    "future-radiance-download-snow-history": PROJECT_ROOT / "src" / "fetch_mod10cm_history.py",
    "future-radiance-monthly-snow-amplification": PROJECT_ROOT / "src" / "test_fixed_nuts2016_monthly_snow_amplification.py",
    "future-radiance-monthly-snow-control-years": PROJECT_ROOT / "src" / "test_fixed_nuts2016_monthly_snow_control_years.py",
    "future-radiance-snow-2021-test": PROJECT_ROOT / "src" / "test_fixed_nuts2016_snow_viirs_2021.py",
    "future-radiance-historical-snow-control": PROJECT_ROOT / "src" / "test_fixed_nuts2016_historical_snow_control.py",
    "future-radiance-albedo-2021-test": PROJECT_ROOT / "src" / "test_fixed_nuts2016_albedo_viirs_2021.py",
    "future-radiance-2021-surface-correction": PROJECT_ROOT / "src" / "build_fixed_nuts2016_2021_surface_correction.py",
    "future-radiance-corrected-minimum-model": PROJECT_ROOT / "src" / "select_fixed_nuts2016_corrected_temporal_model.py",
    "future-radiance-corrected-endpoint-change": PROJECT_ROOT / "src" / "plot_fixed_nuts2016_corrected_endpoint_change.py",
    "habitability-fetch": PROJECT_ROOT / "src" / "fetch_habitability_inputs.py",
    "habitability-build": PROJECT_ROOT / "src" / "build_habitability_metric.py",
    "habitability-historical-validation": PROJECT_ROOT / "src" / "validate_habitability_historical_change.py",
    "habitability-independent-factor": PROJECT_ROOT / "src" / "test_habitability_independent_factor.py",
    "global-calibration-rotation-test": PROJECT_ROOT / "src" / "test_global_calibration_rotation.py",
    "annual-vnl-2017-pixel-test": PROJECT_ROOT / "src" / "test_annual_vnl_2017_pixel_discontinuity.py",
    "annual-vnl-2017-correction-pass2": PROJECT_ROOT / "src" / "test_annual_vnl_2017_correction_pass2.py",
    "annual-vnl-2017-correction-pass3": PROJECT_ROOT / "src" / "test_annual_vnl_2017_correction_pass3.py",
    "paired-density-h-trends": PROJECT_ROOT / "src" / "plot_paired_density_h_trends.py",
    "calibration-offset-test": PROJECT_ROOT / "src" / "test_calibration_offset.py",
    "calibration-country-test": PROJECT_ROOT / "src" / "test_country_calibration_correction.py",
    "calibration-sensitivity": PROJECT_ROOT / "src" / "calibration_sensitivity.py",
}

COMPOSITE_STAGES = {
    "future-radiance-historical": [
        "density-temporal-change",
        "density-demographic-surface",
        "future-radiance-response",
    ],
}


def run_script(script: Path, config: str) -> None:
    command = [sys.executable, str(script), "--config", config]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=[*STAGES, *COMPOSITE_STAGES, "all"])
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    if args.stage in COMPOSITE_STAGES:
        order = COMPOSITE_STAGES[args.stage]
    elif args.stage == "all":
        order = [
            "nuts",
            "validate",
            "radiance",
            "calibration",
            "population",
            "population-recover",
            "merge",
            "fit",
            "trends",
        ]
    else:
        order = [args.stage]

    for stage in order:
        run_script(STAGES[stage], args.config)


if __name__ == "__main__":
    main()
