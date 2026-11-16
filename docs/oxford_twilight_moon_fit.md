# Oxford museum: exploratory solar–lunar DoLP fit

The museum capture on 14 May 2019 at 21:09:51 UTC had Moon elevation/azimuth 41.20/187.15 degrees and Sun elevation/azimuth -10.23/319.08 degrees. Twilight must be considered before interpreting disagreement with a lunar-only pattern as a coordinate error.

Barta et al. (2014), Applied Optics 53, 5193–5204, doi:10.1364/AO.53.005193 measured transitions between solar and lunar patterns and simulated their Stokes addition. Multiple scattering modifies twilight polarization; see Ugolnikov and Maslov, arXiv:1503.01635. These motivate a combined-source test, not a claim that single scattering is accurate at this solar depression.

The script fits the checksum-locked linear RAW DoLP at fixed astronomical source positions and fixed camera registration. Each component has Rayleigh-shaped intensity (1+cos² scattering angle), normalized to unit zenith intensity, and linear polarization proportional to pmax*sin² scattering angle with ideal rear-analyzer transport. This is a geometric proxy. It is not the existing KS91/spectral model and is not calibrated twilight radiative transfer.

The combined I,Q,U are summed before calculating DoLP. The fit estimates a nonnegative solar/lunar zenith-intensity ratio and uniform unpolarized background in lunar zenith-intensity units. Lunar-only and solar-only controls fit their own dilution. The maximum component polarization defaults to a fixed 0.8; it is not fitted, and is not the production atmospheric-factor parameter despite its numerical value. No astronomical position, image rotation, reflection, polarizer zero, per-exposure gain or existing ALAN setting is fitted or changed.

Squared DoLP errors are minimized on two thirds of fixed 96-pixel spatial blocks; the remaining blocks assess prediction. AoP spatial changes are reported independently without calibrating its unknown constant zero. Default fitting uses the central 60-degree full field, rejecting nonphysical measured DoLP and source singularities. An explicit original-shape boolean NPY clear-sky mask can be supplied with --sky-mask. Buildings/clouds are not automatically identified. Maps outside the fitted mask are extrapolations, not validation.

```bash
git pull origin publication-figures

python src/fit_spider_twilight_moon.py \
  --reference-lock data/processed/spider_registration_museum_linear_01/locked_site_references.json \
  --site museum \
  --output data/processed/museum_twilight_moon_fit_01
```

Outputs: twilight_moon_fit.png, twilight_moon_fit.json and twilight_moon_fit_arrays.npz. Inspect held-out errors, parameter bounds, conditioning and the unfitted AoP-change check. A fitted mixture can absorb instrument errors; improvement alone does not identify the physical twilight fraction. Repeating with --pmax 0.5 and 1.0 in separate output directories provides a simple sensitivity check.

This is a capture-specific empirical test. A global solar-elevation switch or general twilight-transition model is deliberately deferred. Production PRC, 1-degree heading calculations and ALAN=5 figures remain unchanged until a model is justified.
