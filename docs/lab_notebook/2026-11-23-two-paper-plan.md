# 2026-11-23 — Two-paper research plan: spectral effects on polarization cues and historical VIIRS interpretation

## Decision and scope

The investigator has decided to develop two connected papers. Paper 1 establishes the biological modeling framework and compares lighting spectra against experimental sky polarization. Paper 2 returns to the corrected historical VIIRS record and asks which biological consequences can be inferred when lighting spectra and operations change.

This entry records an analysis and figure plan, not manuscript prose. The five-figure outline for Paper 1 is the investigator's intended structure. The five-figure outline for Paper 2 below is a provisional organization of the requested analyses. Novelty, intervention attribution and generality of the pilot percentages remain to be established.

**All sites are in scope.** WClawn and Oxford Business Park are the completed spectral-equivalence pilot, not the final study population. Maintain a full-site ledger of RAW availability, metadata, registration, Moon identity/position, orientation, analyzer convention, fitted results and unresolved requirements. Do not infer successful checks or silently exclude pending sites.

## Paper 1 — Modeling the biological consequences of solar–lunar polarization and artificial-light spectra

### Main question

How do the natural solar and lunar contributions, together with HPS or LED skyglow, determine the polarization pattern and receptor contrast available to an animal?

Use the experimental data to constrain and assess the sky pattern, including AoP structure, DoLP-band orientation and PRC. Compare added-light amounts on explicitly stated radiometric and human-photometric bases.

### Five planned figures

| Figure | Purpose | Planned analysis and display |
| --- | --- | --- |
| 1. Natural sky pattern | Establish when lunar-only versus combined solar–lunar modeling is needed | Independently anchored experimental and calculated DoLP/AoP patterns across representative solar/lunar geometries, including twilight and lunar-dominated captures. Show Moon/orientation evidence and spatial agreement. |
| 2. Artificial-light dilution | Reproduce experimental patterns across the site brightness range | Natural patterns plus unpolarized HPS/LED scenarios; experimental/model maps, PRC and receptor-orientation responses. Show how the usable band changes. Apply one consistent workflow across all eligible sites and report pending sites. |
| 3. Equal biological effect | Compare the light amounts required for the same PRC | Fit HPS and LED amounts to the same animal-PRC target, using the agreed green-camera approximation where applicable. Show PRC-versus-band-energy curves, matched-target intersections, percentage differences and sensitivity to spectra, pigment and atmospheric assumptions. |
| 4. Equal human luminance | Compare biological consequences under a human-centered equality criterion | Normalize spectral sky radiances to the same photopic luminance using V(lambda), then calculate animal-weighted PRC. Keep this distinct from equal radiant energy and from equal road illuminance. |
| 5. Documented HPS-to-LED conversion | Apply the model to a real lighting change | Select a confirmed conversion without operational dimming; establish fixture spectra/output, distributions, schedules, footprint and dates where possible. Model the resulting PRC difference and distinguish predicted effects from measured biological change. |

Equal biological effect asks how much light is needed. Equal human luminance asks what PRC results at the same human-weighted brightness.

For Figure 4, equal sky luminance is directly defined from the spectral sky model. Equal road illumination requires additional fixture, geometry and atmospheric-transfer information; do not treat the two as interchangeable.

For Figure 5, "no dimming" is not sufficient evidence for equal emitted flux or unchanged upward emission. Audit lumen output, fixture optics, shielding, installation extent and operating schedules before selecting the case.

### Modeling decisions and work still needed

- Keep astronomical Sun/Moon positions fixed and image registration independently constrained by landmarks and Moon observations. Do not fit a sky rotation or reflection merely to improve pattern agreement.
- Add solar and lunar I, Q and U in a common basis; include the declared artificial-light component; integrate spectral Stokes through the receptor response before calculating polarization/receptor contrast.
- Fit HPS and LED amounts to **PRC**, not to camera DoLP as a substitute for the animal endpoint.
- The investigator accepts camera-green reconstructed PRC as an approximate RH1 target. Retain that explicit approximation and assess it using the published D850 response rather than demanding a unique spectral reconstruction before proceeding.
- Use the ELF D850/Sigma 8 mm calibration from Nilsson/Smolka and colleagues for the camera-response sensitivity comparison. It is not a calibration of this project's added polarizer.
- Separate fitting the natural solar–lunar mixture from comparing HPS/LED amounts. Current pilot source ratios were fitted to DoLP and checked using spatial AoP changes; a joint DoLP/AoP constraint has been proposed but is not yet implemented.
- Current spectral-equivalence calculations use the same transferred natural spectral template for Sun and Moon. A separate twilight spectrum and more complete atmospheric treatment remain sensitivity/development tasks.
- Do not pursue additional PRC-specific adjustments just to remove every residual. Keep the focus on reproducing the pattern and comparing biologically weighted spectral consequences; report residuals and limitations.
- A universal solar-elevation transition rule is not established. Detailed transition modeling may be reserved for a later paper.
- Test spectral alternatives and atmospheric assumptions before promoting one pilot percentage to a general conclusion.
- Check the literature before claiming novelty. The investigator's proposed novelty is a quantitative HPS/LED comparison for a real visual navigational cue, not simply a difference in spectral overlap.

## Pilot evidence at this planning point

These are conditional spectral-model results, not measurements of actual site lamp mixtures.

| Site | Reconstructed green-proxy target maximum absolute PRC | Required HPS amount | Required LED amount |
| --- | ---: | ---: | ---: |
| WClawn | 0.114599 | 6.837084 | 4.166440 |
| Oxford Business Park | 0.111506 | 7.355343 | 4.482261 |

Amounts are uniform added sky radiant energy over **352–700 nm**, relative to lunar zenith radiant energy over the same band. They are not electrical power, lux or VIIRS radiance.

For the current 510 nm RH1 proxy and selected spectral scenarios, LED/HPS required energy = **0.609388**, or **39.061% less LED energy** at the same PRC. Repeating the forward response calculation with the measured D850-green sensitivity and the same numerical target gives **10.459% less LED energy**. This is a response-sensitivity comparison, not a unique correction of the measured target to animal spectral Stokes.

Under the present fixed-natural-field, uniform-unpolarized-addition assumptions, the required LED/HPS ratio is set by the lamp/receptor spectral overlaps. Its equality across the two sites is by construction, not two independent empirical estimates of a universal effect. The site geometry and PRC target determine the absolute amounts required.

At matched RH1 PRC, the HPS and LED RH1 DoLP maps coincide within numerical precision. Display the shared prediction once rather than implying different spatial maps at this matched endpoint. Equal-energy or equal-luminance comparisons answer different questions.

The PRC-matched spatial model has central-field DoLP RMSE **0.02196** at WClawn and **0.02342** at Oxford Business Park, with mean model-minus-measured biases **-0.00867** and **-0.01521**. Matching integrated PRC does not force local DoLP agreement. These are comparisons with the accepted camera-green proxy.

## Paper 2 — Historical VIIRS trends, corrections and biological interpretability

### Main question

What historical changes in artificial light can be established from VIIRS, and how does spectral mismatch limit inference about biologically relevant polarization cues?

Retain the focus on the corrected 2014–2024 historical analysis, including the baseline discontinuity and 2021 snow anomaly. Reconnect the documented LED conversions and dimming examples to the evidence about streetlighting's contribution to total observed emission.

### Provisional five-figure organization

| Figure | Purpose | Planned analysis and display |
| --- | --- | --- |
| 1. Comparable historical record | Establish the effect of necessary corrections | Raw versus baseline-adjusted and snow-corrected trajectories/maps; separate correction effects from substantive lighting change and retain uncertainty/provenance. |
| 2. Sources of observed change | Interpret LED conversions and dimming examples | Evidence-based local comparisons, using intervention timing/footprints and the streetlighting-contribution literature. Do not assume satellite emission is all streetlighting. |
| 3. Spectral ambiguity | Demonstrate why VIIRS alone cannot uniquely recover PRC history | Compare instrument and animal spectral weighting; show HPS- and LED-conditioned biological trajectories and how a spectral transition can alter their interpretation. |
| 4. Nonlinear dimming benefit | Show why starting brightness matters | Apply comparable proportional reductions at different starting radiance/light levels and show the resulting PRC changes under HPS/LED scenarios. Identify the range where cue improvement is largest. |
| 5. Two documented dimming-only interventions | Demonstrate contrasting real-world consequences | Choose two cases in different starting brightness ranges, with independently documented dimming and no concurrent spectral conversion in the relevant comparison. Show observed radiance changes and conditional PRC consequences. |

### Interpretation to preserve

- VIIRS radiance alone does not uniquely identify a PRC trajectory when source spectra change. A decrease in measured radiance need not be a uniquely determined improvement in the animal cue.
- Stable spectra, documented operational changes and agreement across alternative spectral scenarios can still support useful directional conclusions.
- Show HPS-only and LED-only PRC scenario envelopes. Treat them as conditional scenarios, not guaranteed extrema over every possible spectrum, atmosphere and spatial lighting distribution.
- The two dimming cases should illustrate the nonlinear response: the same proportional light reduction can produce different PRC improvements at different starting levels.
- Match intervention footprints, observation periods and aggregation scales; document other emitting sectors and concurrent changes.
- The broader implication is the need for biologically informative spectral measurements and instrument response functions, rather than a direct universal radiance-to-PRC conversion.
- The previously discussed paper about how much streetlights contribute is an evidence strand to identify and cite accurately, not an assumed universal streetlighting fraction.

## Order of work and deliverables

1. **Complete the all-site data/check ledger.** Preserve raw files and exclusions, resolve pending geometry/Moon checks, and record why any site cannot yet be used.
2. **Settle the natural-pattern model and comparison criteria.** Keep fixed geometry; assess both DoLP and offset-independent AoP structure. State whether a PRC endpoint is maximum, heading median or a specified receptor orientation; do not mix endpoints.
3. **Run the same HPS/LED PRC-equivalence analysis across eligible sites.** Keep identical supports, explicit band-energy normalization, source hashes and pending-site accounting.
4. **Expand spectral and atmospheric sensitivity analyses.** Evaluate representative LED spectra, HPS, receptor sensitivity and the treatment of twilight. Assess robustness rather than optimizing away remaining residuals.
5. **Implement the equal-photopic-luminance comparison.** Use a documented V(lambda) dataset, consistent units and normalization; keep energy, luminance and road illuminance separate.
6. **Audit intervention candidates.** Select one conversion-without-dimming example and two dimming-only examples based on independent documentary evidence before interpreting satellite changes.
7. **Return to the corrected historical analysis.** Use the established model to express biological consequences conditionally and compare dimming in different parts of the nonlinear PRC response.
8. **Finalize figures after analysis choices are stable.** Main figures should communicate the biological comparisons; retain residuals, calibration checks and full site results as companion/supplementary material.

No additional production analyses or historical corrections are executed by this planning entry.

## Current implementation and restart pointers

- [PRC-equivalent lamp scenarios and D850-response provenance](../prc_lamp_equivalence.md)
- [Camera-band PRC comparison](../mixed_prc_comparison.md)
- [Intensity-weighting diagnostic](../prc_weighting_audit.md)
- [Figure design notes](../prc_equivalence_figure_design.md)
- [Signed analyzer convention](../signed_analyzer_convention.md)
- [Existing historical pipeline notes](../pipeline.md)

Relevant code:

- `src/fit_spider_twilight_moon.py`
- `src/audit_spider_mixed_aop.py`
- `src/compare_spider_mixed_prc.py`
- `src/fit_spider_prc_lamp_equivalence.py`
- `src/plot_spider_prc_equivalence_figure.py`
- `src/plot_spider_prc_equivalence_batch.py`

EcoVis pilot inputs/results:

- `data/processed/spider_sign_checked_01/registration/locked_site_references.json`
- `data/processed/spider_sign_checked_01/<site>/fit`
- `data/processed/spider_prc_checked_01/<site>`
- `data/processed/spider_uv510_01/inputs` — the selected spectral-input set; do not search all processed versions indiscriminately
- `data/processed/spider_prc_lamp_equivalence_01/<site>`
- `data/processed/spider_prc_publication_01`

The two-site pilot figures are layout and analysis prototypes. Batch plotting supports completed results; it does not itself complete missing RAW reconstruction, registration, Moon checks or spectral fitting for other sites.
