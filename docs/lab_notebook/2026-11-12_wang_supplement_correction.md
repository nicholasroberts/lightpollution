# Wang supplement completion and coordinate checks — 2026-11-12

Read investigator-supplied jopt065601-suppdata.pdf, SHA256 489b3738e4cf48a6841dc2daff692ee94a5feeefe860631a27479e729f9e7d93; visually inspected page 4. S4/S5 duplicate the y-type matrix. Investigator accepted typo explanation and authorized using distinct original Preetham matrices. Implemented corrected x/y using cited Preetham author manuscript Appendix A.2, while preserving supplement_as_printed as an explicitly selectable branch. No silent change to printed source. Sources: https://thinks.github.io/docs/1999_a_practical_analytic_model_for_daylight.pdf (Appendix A.2, pp22-23, Table2 p28), DOI 10.1145/311535.311545. Coefficients reflect this author-manuscript version, not an assertion that every published variant has identical rounding.

Added Wang supplemental S1-S3 distribution matrices; retained D_Y intercept -2.5571 exactly as supplied (Preetham manuscript has -2.5771), not silently corrected. Added 41-row 380..780nm S0/S1/S2 basis from Table2. Numerical interpolation inside covered range only; UV below380 and NIR beyond780 not fabricated. This coverage still cannot support a complete future VIIRS/HPS conversion on its own.

Implemented zenith chromaticity, normalized directional x/y (eq10/11), relative daylight spectra (eq13-15) and eq16 factor. The 90deg-reference ray is explicitly in source-zenith meridian opposite source azimuth, theta90=pi/2-theta_source. This is a stated reference choice because a scattering angle alone does not uniquely select a sky ray. Ssource uses circumsource skylight chromaticity, not a validated direct solar/lunar disc spectrum. Moon spectral substitution remains a separate physical task.

Full printed eq17 audit at source elevation45deg, azimuth30deg, T3, assumed neutral separation24deg, wavelengths380/450/550/780: corrected zenith x=.245678338,y=.251476440; printed duplicate x=y=.252170038. No nonphysical DoLP pixels in these eight cases. Corrected maximum DoLP approx .2014/.1999/.1996/.2172. This is a synthetic numerical audit, not reproduction of paper figures or empirical skies. Raw eq16 S values include negatives under the declared spectral/reference construction; not claimed to be a universal 0..1 factor. No output clipped to force agreement.

Coordinate audit: simultaneous source/view azimuth rotation should preserve DoLP and transform a common global-plane Q+iU by spin2. Literal printed eq4/5 plus eq3 has average DoLP covariance residual .0005003 and complex spin2 residual .208983 in the specified synthetic sample. This may reflect reference-basis conventions or transcription/parameterization issues; do NOT interpret as a proven paper error. Added source_meridian_pattern as an explicitly separate comparison: evaluate source on +y with eq2/3 and rotate common complex field by exp(2i angle). Covariance residual ~5.5e-17. This convention is NOT silently substituted into default natural_stokes or earlier pilot. Needs actual image-basis validation.

14 model tests pass, including corrected x/y separation, zenith normalization,560nm basis anchor, spectral source anchor, no extrapolation and explicit meridian covariance. New audit command produces spectral_coefficient_audit.png, aop_coordinate_audit.png, wang_spectral_audit.json. Charts generated and inspected locally. No actual lunar ephemeris/photograph fit run. Existing dynamic pilot stays explicitly reduced until convention/reference choices verified.

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/audit_wang2016_spectral_model.py \
  --neutral-separation-deg 24 \
  --output data/processed/wang_spectral_coordinate_audit_01
```

Upload both PNGs and JSON. Keep illustrative24deg labelled as assumed. Next compare AoP reference-frame maps to observed skies, verify capture UTC/orientation and only then couple full spectral model to true Moon tracks. Locked PRC data/fits and prior results untouched.
