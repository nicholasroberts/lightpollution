#!/usr/bin/env python3
"""Calculate bounded all-LED spectral endpoints for future VIIRS predictions.

This stage does NOT forecast a universal country intercept. Instead it isolates
one physically bounded component of the observed radiance level: the change in
VIIRS response caused by a shift in emitted spectrum at fixed green-band
radiance.

Two published relationships are used:

1. Sánchez de Miguel et al. (2019), Remote Sensing of Environment 224,
   Table 2, empirical DSLR synthetic-photometry relation between correlated
   colour temperature (CCT) and Nikon-D3s-like ISS G/R:

      CCT / 1e4 =
          -3.0 x^4 + 5.8 x^3 - 3.2 x^2 + 1.0 x + 0.06
      where x = G/R and the published valid interval is 0.2 <= G/R <= 1.

2. Sánchez de Miguel et al. (2022), Science Advances 8:eabl6891,
   their Eq. 2:

      G_ISS / VIIRS = 0.21 * 1.5 ** (G_ISS / R_ISS)

For equal G-band output, VIIRS/G is therefore the reciprocal of Eq. 2.
The ratio of future to present VIIRS response at equal G output is:

      VIIRS_future / VIIRS_present
          = (G/VIIRS)_present / (G/VIIRS)_future

The current spectral anchor is the post-period country G/R median from the
processed ISS analysis (2014-2020 composite), not the 2024 regression
intercept. Country intercept histories remain separate and can continue to
contain non-spectral effects.

The endpoint scenarios are deliberately bounded all-LED CCTs (default
2700/3000/4000/5000 K). These are scenario endpoints, not claims that every
street or building light will converge to one CCT.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from project_config import configured_path, load_config


ARTICLE_2019_DOI = "10.1016/j.rse.2019.01.035"
ARTICLE_2022_DOI = "10.1126/sciadv.abl6891"


def cct_from_gr(gr: float | np.ndarray) -> float | np.ndarray:
    """Published CCT relation; returns kelvin."""
    x = np.asarray(gr, dtype=float)
    return 1.0e4 * (
        -3.0 * x**4
        + 5.8 * x**3
        - 3.2 * x**2
        + 1.0 * x
        + 0.06
    )


def gr_from_cct(cct_k: float, lo: float = 0.2, hi: float = 1.0) -> float:
    """Invert the published polynomial on its stated valid G/R interval."""
    target = float(cct_k)
    flo = float(cct_from_gr(lo) - target)
    fhi = float(cct_from_gr(hi) - target)
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if flo * fhi > 0:
        raise ValueError(
            f"CCT {target:g} K lies outside the invertible published "
            f"G/R interval [{lo}, {hi}]."
        )
    return float(brentq(lambda x: cct_from_gr(x) - target, lo, hi))


def g_over_viirs_from_gr(gr: float | np.ndarray) -> float | np.ndarray:
    """Sánchez de Miguel et al. (2022), Eq. 2 central relation."""
    x = np.asarray(gr, dtype=float)
    return 0.21 * np.power(1.5, x)


def viirs_per_g_from_gr(gr: float | np.ndarray) -> float | np.ndarray:
    return 1.0 / g_over_viirs_from_gr(gr)


def build_led_scenarios(cfg: dict) -> pd.DataFrame:
    ccts = [float(x) for x in cfg["led_cct_scenarios_k"]]
    labels = cfg.get("led_cct_labels", {})

    rows = []
    for cct in ccts:
        gr = gr_from_cct(cct)
        g_over_viirs = float(g_over_viirs_from_gr(gr))
        viirs_per_g = 1.0 / g_over_viirs
        rows.append(
            {
                "scenario": labels.get(str(int(cct)), f"{int(cct)} K LED"),
                "cct_k": cct,
                "endpoint_gr": gr,
                "endpoint_g_over_viirs": g_over_viirs,
                "endpoint_viirs_per_g": viirs_per_g,
            }
        )

    table = pd.DataFrame(rows).sort_values("cct_k").reset_index(drop=True)

    baseline_cct = float(cfg["legacy_reference_cct_k"])
    baseline_gr = gr_from_cct(baseline_cct)
    baseline_vg = float(viirs_per_g_from_gr(baseline_gr))

    table["legacy_reference_cct_k"] = baseline_cct
    table["legacy_reference_gr"] = baseline_gr
    table["viirs_multiplier_vs_legacy_equal_g"] = (
        table["endpoint_viirs_per_g"] / baseline_vg
    )
    table["pct_viirs_change_vs_legacy_equal_g"] = 100.0 * (
        table["viirs_multiplier_vs_legacy_equal_g"] - 1.0
    )
    table["delta_log10_viirs_vs_legacy_equal_g"] = np.log10(
        table["viirs_multiplier_vs_legacy_equal_g"]
    )
    return table


def build_country_plateaus(
    config: dict,
    scenarios: pd.DataFrame,
) -> pd.DataFrame:
    processed = configured_path(config, "processed")
    source = processed / "iss_country_spectral_shift.csv"
    if not source.exists():
        raise FileNotFoundError(
            f"Missing {source}. Run: "
            "python src/run_pipeline.py iss-spectral-validation"
        )

    country = pd.read_csv(source)
    if "post_gr_median" not in country.columns:
        raise KeyError(f"{source} has no post_gr_median column.")

    country = country[
        ["CNTR_CODE", "post_gr_median", "n_common_valid_pixels"]
    ].copy()
    country = country.rename(
        columns={
            "post_gr_median": "spectral_anchor_gr",
            "n_common_valid_pixels": "spectral_anchor_common_pixels",
        }
    )
    country["spectral_anchor_implied_cct_k"] = cct_from_gr(
        country["spectral_anchor_gr"]
    )
    country["spectral_anchor_g_over_viirs"] = g_over_viirs_from_gr(
        country["spectral_anchor_gr"]
    )
    country["spectral_anchor_viirs_per_g"] = viirs_per_g_from_gr(
        country["spectral_anchor_gr"]
    )

    rows = []
    for crow in country.itertuples(index=False):
        if not np.isfinite(crow.spectral_anchor_gr):
            continue

        for srow in scenarios.itertuples(index=False):
            multiplier = (
                crow.spectral_anchor_g_over_viirs
                / srow.endpoint_g_over_viirs
            )
            rows.append(
                {
                    "CNTR_CODE": crow.CNTR_CODE,
                    "scenario": srow.scenario,
                    "endpoint_cct_k": srow.cct_k,
                    "spectral_anchor_gr": crow.spectral_anchor_gr,
                    "spectral_anchor_implied_cct_k": (
                        crow.spectral_anchor_implied_cct_k
                    ),
                    "endpoint_gr": srow.endpoint_gr,
                    "remaining_delta_gr_to_endpoint": (
                        srow.endpoint_gr - crow.spectral_anchor_gr
                    ),
                    "spectral_anchor_common_pixels": (
                        crow.spectral_anchor_common_pixels
                    ),
                    "spectral_anchor_g_over_viirs": (
                        crow.spectral_anchor_g_over_viirs
                    ),
                    "endpoint_g_over_viirs": srow.endpoint_g_over_viirs,
                    "viirs_multiplier_at_equal_g": multiplier,
                    "pct_viirs_change_at_equal_g": 100.0 * (multiplier - 1.0),
                    "delta_log10_viirs_at_equal_g": float(
                        np.log10(multiplier)
                    ),
                }
            )

    return pd.DataFrame(rows)


def make_plots(
    scenarios: pd.DataFrame,
    country: pd.DataFrame,
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.plot(
        scenarios["cct_k"],
        scenarios["pct_viirs_change_vs_legacy_equal_g"],
        marker="o",
    )
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("All-LED endpoint CCT (K)")
    ax.set_ylabel("VIIRS response change vs legacy reference (%)")
    ax.set_title(
        "Bounded spectral effect on VIIRS at equal green-band output"
    )
    fig.tight_layout()
    out = figure_dir / "led_endpoint_viirs_vs_cct.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Wrote: {out}")

    # Country-by-scenario remaining correction.
    pivot = country.pivot(
        index="CNTR_CODE",
        columns="endpoint_cct_k",
        values="pct_viirs_change_at_equal_g",
    )
    pivot = pivot.sort_index()

    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    x = np.arange(len(pivot))
    for cct in pivot.columns:
        ax.plot(
            x,
            pivot[cct],
            marker="o",
            linewidth=1.0,
            label=f"{int(cct)} K",
        )
    ax.axhline(0, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=90)
    ax.set_ylabel(
        "Remaining VIIRS spectral correction from ISS post-period (%)"
    )
    ax.set_xlabel("Country")
    ax.set_title(
        "Country-specific bounded all-LED spectral endpoint\n"
        "equal green-band output; not a forecast of total intercept change"
    )
    ax.legend(title="LED endpoint")
    fig.tight_layout()
    out = figure_dir / "led_endpoint_country_remaining_correction.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Wrote: {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/pipeline.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    cfg = config["analysis"]["led_spectral_plateau"]

    scenarios = build_led_scenarios(cfg)
    countries = build_country_plateaus(config, scenarios)

    processed = configured_path(config, "processed")
    scenario_path = processed / "led_spectral_plateau_scenarios.csv"
    country_path = processed / "led_country_spectral_plateaus.csv"
    summary_path = processed / "led_spectral_plateau_summary.txt"

    scenarios.to_csv(scenario_path, index=False)
    countries.to_csv(country_path, index=False)

    figure_dir = (
        configured_path(config, "figures")
        / "population_relationship"
        / "spectral_validation"
    )
    make_plots(scenarios, countries, figure_dir)

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("BOUNDED ALL-LED SPECTRAL ENDPOINTS FOR VIIRS\n")
        handle.write("=" * 72 + "\n\n")
        handle.write(
            f"CCT->G/R source: Sánchez de Miguel et al. 2019, "
            f"DOI {ARTICLE_2019_DOI}\n"
        )
        handle.write(
            f"G/R->G/VIIRS source: Sánchez de Miguel et al. 2022, "
            f"DOI {ARTICLE_2022_DOI}\n"
        )
        handle.write(
            "Interpretation: all VIIRS multipliers are spectral-only changes "
            "at equal ISS green-band output. They are not forecasts of total "
            "night-time radiance or of the full country intercept.\n\n"
        )
        handle.write("ALL-LED ENDPOINT SCENARIOS\n")
        handle.write(
            scenarios.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        handle.write("\n\nCOUNTRY-SPECIFIC REMAINING SPECTRAL CORRECTIONS\n")
        handle.write(
            countries.to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}",
            )
        )
        handle.write("\n\nMODEL USE\n")
        handle.write(
            "For a future country intercept model, keep the country-specific "
            "non-spectral level/trend separate. Add only a bounded spectral "
            "term that moves from zero at the chosen anchor toward "
            "delta_log10_viirs_at_equal_g. Do not force countries to share "
            "one intercept trajectory.\n"
        )
        handle.write(
            "A later forecasting stage can parameterize the approach to the "
            "plateau with a bounded adoption fraction f_LED(t) in [0,1].\n"
        )

    print("\nALL-LED ENDPOINT SCENARIOS")
    print(
        scenarios[
            [
                "scenario",
                "cct_k",
                "endpoint_gr",
                "endpoint_g_over_viirs",
                "pct_viirs_change_vs_legacy_equal_g",
                "delta_log10_viirs_vs_legacy_equal_g",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nCOUNTRY-SPECIFIC REMAINING CORRECTION")
    print(
        countries[
            [
                "CNTR_CODE",
                "scenario",
                "spectral_anchor_gr",
                "spectral_anchor_implied_cct_k",
                "endpoint_gr",
                "pct_viirs_change_at_equal_g",
                "delta_log10_viirs_at_equal_g",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6g}",
        )
    )

    print("\nWROTE")
    for path in [scenario_path, country_path, summary_path]:
        print(path)


if __name__ == "__main__":
    main()
