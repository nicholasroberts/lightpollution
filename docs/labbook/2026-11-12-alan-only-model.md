# ALAN-only depolarization — 12 November 2026

Investigator decision: atmospheric turbidity is not an input to the biological model. Added unpolarized ALAN increases I and leaves Q/U unchanged. For a spatially constant relative addition a=I_ALAN/I_natural, p_total=p_natural/(1+a). An equivalent dimensionless ALAN depth t=ln(1+a) gives exp(-t); it is not atmospheric turbidity and has no Wang eq7 intercept.

The active measured comparison now uses the reduced Wang angular/horizon structure (eq6+eq8) with no eq7 multiplier. The previous full eq17 spectral extension is omitted because its Perez/Preetham spectra themselves depend on atmospheric turbidity. This is a declared model change, not an exact reproduction of the previous spectral model with a single factor removed. The spectral and atmospheric functions remain available for archival equation audits, not active biological inference. The pilot also defaults to no atmospheric multiplier.

No-ALAN does not establish the theoretical maximum possible DoLP: the retained angular structure, horizon term and neutral separation still constrain its amplitude. In particular, the prior estimate of Cottage mean 0.334 was the old spectral result divided by eq7 attenuation; it is not a prediction of this new reduced model. Do not interpret a fitted ALAN amount until the reduced baseline is checked against measurements. No Cottage amplitude calibration has been applied.

The relative-addition option assumes ALAN proportional to natural intensity spatially. A uniform absolute skyglow instead needs the existing array-valued add_unpolarized operation and generally gives position-dependent dilution. Neither option converts VIIRS radiance into ALAN without further calibration.

Saved alignments and Moon anchors are reused. The 60-degree full field remains 30-degree radius. AoP still has one fitted constant coordinate offset. Historical PRC is not a fitting target.

Validation: 22 Wang-related numerical tests passed, including no-atmospheric default, unchanged Q/U, exact exp(-t) dilution and rejection of negative ALAN depth. Actual Cottage/Oxford measurements require execution on EcoVis; no new measured result is claimed.
