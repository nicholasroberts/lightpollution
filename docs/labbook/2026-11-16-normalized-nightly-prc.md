# Normalized working model: spectral sky maps and nightly spider PRC

The investigator selected the Rayleigh-denominator candidate after the matched
Cottage/Skillinge/Slabhouse shape comparison. It now becomes the working default
for the spider spectral/nightly forward workflow. Its finite-neutral-splitting
extension remains a hypothesis, not a full multiple-scattering solution or the
complete Wang eq17. Keep the previous envelope as a comparison option.

## Common model across maps and receptors

New spider_natural_polarization.py supplies complex fractional polarization:

P_natural = W_Berry * H(theta)/(1+cos²gamma)  [rayleigh_normalized]
P_natural = W_Berry * H(theta)                [berry_envelope, previous law]

The phase and neutral points are identical between laws. The helper includes
the horizon factor ONCE. Callers apply the residual atmospheric multiplier eta
ONCE to Q/U. KS91 provides natural relative intensity; unpolarized ALAN adds to
I only. The resulting Stokes signals are integrated spectrally and spatially
into orthogonal receptor catches before their log ratio is taken.

The helper is used by nightly field integration, capture-time spectral maps,
and the optional historical Cottage scenario-adjustment script. Archival
wang2016_sky.py equation-audit functions are deliberately unchanged. No empirical
Europe radiance–PRC model or raw measurement is changed.

--natural-dop-law rayleigh_normalized is the new nightly/runner default;
--natural-dop-law berry_envelope explicitly selects the old amplitude law.
New nightly/runner default neutral separation is10deg, matching the measured
comparison (older nightly scripts defaulted to24deg). For an isolated amplitude-law
comparison, run both options with the SAME explicit10deg separation. Historical
24deg runs require both --natural-dop-law berry_envelope and
--neutral-separation-deg24. Separation is still a prescribed assumption.

## Freeze the existing ALAN scenario

The runner accepts --alan-configs DIR. It copies HPS_site_alan.json and
LED_mid_site_alan.json byte-for-byte and records hashes. Reuse mode disallows
Cottage retargeting, a new skyglow ratio or regenerated DNB response. It does
not read experimental DoLP to set ALAN. Existing source-scale assumptions,
including any earlier user-requested Cottage0.45 target, remain explicit
provenance. That earlier target is NOT reimposed under the new law.

Those configurations are scenarios, not measured capture-time ground skyglow.
In particular a reused Cottage target adjustment is not newly independent
calibration merely because it is frozen. Titles call it a prior scenario rather
than claiming the current mean equals0.45. Actual new FoV mean/median values
are calculated and annotated on the maps.

Spectra remain the existing352–750nm inputs, green510nm alpha-band sensitivity,
HPS and LED templates, and recorded DNB-overlap scaling. Body tilt stays0deg,
heading0–350deg every10deg, field60deg full width,1024 field samples and20min
night sampling. Eta=.8 is retained as a prescribed sensitivity value; it was
not fitted to the observed reduction. Level0 provides a natural-only control.

## EcoVis command: normalized model, one output folder

Reuse inputs and configurations from the preceding UV/green510 scenario run:

```bash
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/run_independent_site_alan.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --spectra data/processed/spider_uv510_01/inputs \
  --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
  --alan-configs data/processed/spider_uv510_01/results/alan_configs \
  --natural-dop-law rayleigh_normalized \
  --neutral-separation-deg 10 \
  --atmospheric-factor 0.8 --levels 0 5 6 7 \
  --output data/processed/spider_normalized_prc_01
```

The spectral/nightly pilot currently covers Cottage and Oxford Business Park,
the established low/high-ALAN comparison. Skillinge and Slabhouse supplied the
natural-law validation; extending their nightly skyglow scenarios requires
site source-scale assumptions and is not silently inferred from their DoLP.

Look in **data/processed/spider_normalized_prc_01**:

- comparison_level_0/visible_spectral_nightly_comparison.png: natural control,
  sky maps with measured camera panels, local FoV means/medians and heading/time
  summary. With the same natural spectrum and sensitivity, HPS/LED agree at
  level0; differences then arise from the added lamp spectra.
- comparison_level_5/, comparison_level_6/, comparison_level_7/: same layout
  under fixed source additions. LED/HPS spectral contributions remain distinct.
- hps/ and led/: nightly signed PRC heatmaps, absolute-PRC heading ranges,
  heading_responses.csv, prc_ranges.csv and provenance.json.
- scenario_notes.json and alan_configs/: frozen source assumptions and hashes.

Rerun the identical command with --natural-dop-law berry_envelope and a fresh
--output data/processed/spider_previous_law_prc_01 to compare laws without
altering source additions, atmosphere or neutral splitting. No original output
folder is overwritten. Differences in absolute PRC across these scenarios are
model sensitivity, not new behavioral thresholds or calibrated detection odds.

## Reproducibility and compatibility

Nightly provenance records the law, separation and shared helper hash. Snapshot
rendering reads that law from provenance. Old provenance without a law field
renders with berry_envelope, so old results are not relabelled or redrawn with
the new law. HPS/LED plotting rejects different laws or neutral separations.
The old camera-conjugation offset is not applied to the physical sky model:
it was a camera-coordinate comparison convention, not an animal-world basis.

## Checks

Five new tests check shared phase/amplitude, explicit/default field behavior,
invalid-law rejection, legacy snapshot provenance and unchanged ALAN config bytes without retargeting. The four existing nightly
checks cover world-rotation invariance, tangent field transport, spectral
normalization and ALAN/atmosphere effects. Coarse full-night executions of both laws check
CSV/JSON and figure creation; it is a software smoke test, not a new validated
site prediction. Real measured comparison panels require EcoVis saved NPZ files.

The aim remains species-weighted PRC ranges across heading and time under
realistic spectral ALAN scenarios. Natural-law refinement is now carried into
that workflow rather than pursued as a separate sky-width fitting exercise.
