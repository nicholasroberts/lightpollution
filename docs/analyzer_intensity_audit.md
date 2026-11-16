# Analyzer-intensity audit

The production linear reconstruction currently uses I=I0+I90, Q=I0-I90, U=I45-I135 after image registration. It assigns nominal analyzer angles. Museum landmark rotations are approximately 0, 42.8, 86.1 and 132.8 degrees. Registration moves pixel locations; it does not itself update the analyzer design matrix.

This diagnostic reloads checksum-verified linear camera-blue TIFFs and applies the saved registration once. It estimates per-pixel I/Q/U for each available complete quartet using two fixed design matrices: nominal angles and relative landmark roll angles. No source geometry, polarizer zero, gains or offsets are fitted. The landmark-roll case is an approximation, not a replacement for field-dependent rear-analyzer transport through the fisheye.

The forward equation is J(theta) = (I + Q cos(2 theta) + U sin(2 theta))/2. Least squares uses all four measurements equally; unlike the existing production formula, both analyzer pair sums inform I. Each omitted-frame check solves from three frames and predicts the fourth. With four measurements and three parameters these checks share one closure constraint and are correlated; they are not independent evidence of physical calibration. Reversing every angle is intensity-degenerate with reversing U and cannot be resolved by closure alone.

The output reports pair-sum closure, forward residuals, omitted-frame errors and disagreement between the two complete quartets. Central 60-degree field samples are common to every exposure; block means are formed after registration. No saturation calibration or sky/structure segmentation is inferred. Statistics are unclipped; plot limits are display choices. A sky changing during the sequence can also cause residuals.

```bash
git pull origin publication-figures

python src/audit_spider_analyzer_intensities.py \
  --reference-lock data/processed/spider_registration_museum_linear_01/locked_site_references.json \
  --site museum \
  --output data/processed/museum_analyzer_intensity_audit_01
```

Return `analyzer_consistency.png` and `analyzer_consistency.json`. No production reconstruction or theoretical values are changed. Absolute analyzer orientation and optical Mueller calibration remain unresolved even if closure improves.
