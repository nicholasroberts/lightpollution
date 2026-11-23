# Current decisions and open questions

## Adopted

- All sites remain in scope; keep pending data/checks visible.
- Preserve legacy scripts and source products; new workflow is additive.
- Astronomical Sun/Moon coordinates and image registration are independently constrained.
- Analyzer sequence sign, camera parity, analyzer zero and sky-view display orientation are distinct.
- Integrate spectral and spatial I/Q/U before forming receptor contrast.
- Fit lamp amount to PRC; use camera green as the accepted approximate RH1 target.
- Primary pilot endpoint is maximum absolute PRC in the zenith-centred 60-degree full field.
- HPS/LED amount is added sky band radiant energy, not electrical power or lux.
- Uniform unpolarized pollution is the current simplifying model.
- Publication sky view: north up, east left; keep numerical basis transformations explicit.
- Main comparison: measurement, shared PRC-matched model and energy curves; residuals separate.

## Pending

- Completion and checks for every site.
- Any joint solar/lunar fitting objective involving spatial AoP.
- Distinct twilight solar spectrum and atmospheric sensitivity treatment.
- Broader HPS/LED spectral scenarios and uncertainty in response assumptions.
- Equal-luminance V(lambda) dataset and its implementation.
- Documentary selection of the HPS-to-LED conversion without dimming.
- Final novelty assessment and figure acceptance.

## Interpretation to preserve

The common pilot LED/HPS energy ratio follows from shared lamp/receptor overlap
under uniform unpolarized addition. Site geometry and target set absolute
amounts. Do not treat equal ratios as independent experimental replication.
D850 green and RH1 sensitivity runs hold the numerical target fixed; they do
not uniquely convert the measured camera target to animal Stokes.

Record later changes in dated lab-book entries with reasons and affected runs.
