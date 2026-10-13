# Investigation record — VIIRS 2017 calibration discontinuity

**Status:** KEEP / production decision reached  
**Decision date:** 2026-10-12  
**Branch:** `future-population-models`

## Scientific question

The Annual VNL VIIRS DNB series shows a known discontinuity around the 2016-to-2017 calibration transition. The project initially carried a correction estimated from a small comparison subset of countries. That was recognized as wrong in scope because a sensor calibration change is a property of the data stream, not of those countries.

The investigation therefore asked:

1. does the discontinuity remain in pixels already retained as lit, or is it created by EOG zero/background masking?
2. is the step more consistent with an additive zero-point shift or a multiplicative change?
3. what correction value and domain are supported well enough for production use?
4. can that correction be applied before NUTS3 aggregation while preserving raw data and exact area weighting?

## Evidence chain

### 1. Retained-pixel test

Stage:

```bash
python src/run_pipeline.py annual-vnl-2017-pixel-test
```

Pixels positive in all of 2015-2018 were followed directly. The diagnostic contrast was:

```text
excess_2017 =
    (R2017 - R2016)
    - 0.5 * [(R2016 - R2015) + (R2018 - R2017)]
```

The positive 2017 excess remained on persistent lit support.

**Conclusion:** EOG background zeroing does not remove the calibration discontinuity. The step is not merely a support/mask artifact.

### 2. Independent-brightness pass

Stage:

```bash
python src/run_pipeline.py annual-vnl-2017-correction-pass2
```

To avoid circularity, brightness class was defined from:

```text
mean(R2013, R2014)
```

rather than from 2016, which also appears in the transition contrast.

The independent-bin excess was approximately flat at about 0.23-0.24 nW cm^-2 sr^-1 through the dim-to-moderate range. The empirical pass-2 value was:

```text
0.23478 nW cm^-2 sr^-1
```

### 3. Segmented 2013-2019 time-series pass

Stage:

```bash
python src/run_pipeline.py annual-vnl-2017-correction-pass3
```

Primary model:

```text
R_t = a + b * (year - 2016) + delta * I(year >= 2017)
```

2013-2019 was used so 2020 could not determine the calibration estimate. Additive and multiplicative formulations were compared. Candidate offsets included literature sensitivities, the pass-2 empirical value and brightness-dependent variants.

The additive form was the more useful representation in the radiance range most sensitive to the zero-point shift. The universal pass-2 value performed well through independent baseline radiance <=10 nW cm^-2 sr^-1.

## Production decision

From 2017 onward:

```text
raw <= 0 : corrected = raw
raw > 0  : corrected = max(raw - 0.235, 0)
```

The rounded **0.235 nW cm^-2 sr^-1** offset is applied:

- to every analysis country/territory;
- to every positive pixel;
- before NUTS3 aggregation;
- with clipping at zero;
- without modifying the original EOG raster.

The production mean then uses the normal exact polygon overlap and spherical physical-area weighting.

## Why this is a pixel correction

A fixed pixel offset does not translate to a fixed NUTS3 offset because regions contain different fractions of positive pixels and different fractions may clip to zero. The pipeline therefore records:

- raw mean;
- corrected mean;
- raw-positive area fraction;
- clipped area fraction;
- effective NUTS3 offset;
- calibration provenance.

Subtracting 0.235 directly from a NUTS3 mean would be incorrect.

## Positive controls and later evidence

The 2021 anomaly investigation supplied a useful later positive control. In the exact mask/intensity decomposition, the known 2017 transition behaved overwhelmingly as an **intensity** shift, whereas the 2021 event contained both support and intensity components. This supports the interpretation that the 2021 event is not simply another copy of the 2017 zero-point change.

## Limitations retained

- The high-brightness end does not follow one perfectly constant additive step.
- The correction is an empirical production harmonization, not a complete physical sensor model.
- Raw EOG values remain preserved for audit.
- Any future revision must reproduce or improve the persistent-pixel, independent-brightness and segmented-time-series diagnostics before replacing the production correction.

## KEEP conclusion

The current production calibration is:

```text
2017+ positive pixels: subtract 0.235 nW cm^-2 sr^-1, clip at zero
```

This result is considered settled enough for production use unless contradictory direct pixel evidence emerges.

## Key commits

```text
e3b5638f  Add Annual VNL 2017 retained-pixel diagnostic
919c2550  Add independent-brightness 2017 calibration pass
6a71efa8  Add segmented time-series VIIRS calibration pass
959e332f  Use Europe-wide pixel-level VIIRS calibration
0642cdf5  Set production VIIRS correction to 0.235 pixel offset
65036193  Promote pixel-level calibration to production pipeline
```

See also `docs/pipeline.md` for implementation details.
