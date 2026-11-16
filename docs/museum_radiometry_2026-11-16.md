# Museum Moon–DoLP reconstruction audit, 16 November 2026

## Finding and scope

Good image-edge registration does not establish a correct polarization reconstruction. The measured broad DoLP band remains differently oriented from the simple lunar-scattering expectation after replacing nominal analyzer angles by measured global image rolls, and after an idealized ray-dependent rear-analyzer transport test. No sky pattern was rotated or fitted to obtain these results. The real optical Mueller response is not calibrated by the transport test.

The supplied 8-bit TIFFs are rendered RGB images. All 12 contain the same Adobe RGB (1998) ICC profile, whose rTRC/gTRC/bTRC tags each specify a one-element `curv` with exponent 2.19921875. XMP records NEF-to-DNG-to-TIFF conversion through Camera Raw and Photoshop, `CameraProfile=Adobe Standard`, and a nested `Look=Adobe Color`. The old `Medium Contrast` curve also appears in metadata, but it is not assumed active: ProcessVersion 11 explicitly records a linear top-level PV2012 curve. The nested Adobe Color look includes a non-identity PV2012 curve.

More directly, the XMP white-balance settings vary between the analyzer triplets:

| Exposures | Nominal analyzer | Temperature (K) | Tint |
|---|---:|---:|---:|
| 1–3 | 0° | 6400 | +6 |
| 4–6 | 45° | 6000, 5900, 5850 | −17, −17, −18 |
| 7–9 | 90° | 6150, 6200, 6250 | −16, −16, −15 |
| 10–12 | 135° | 6400 | +2, +3, +4 |

The blue channel of these rendered exports therefore cannot simply be assumed to have a common linear camera response. Such differences can introduce false Stokes differences. This is a demonstrated processing inconsistency and a plausible contributor to the pattern mismatch, not proof that it is the sole cause.

## Controlled diagnostic tests

All tests use the fixed landmark-derived registration and image Moon position. Astronomical coordinates, ALAN, spectra and theoretical DoLP are unchanged. On the subsampled central 60-degree full-field grid:

| Input interpretation / analyzer model | Median DoLP |
|---|---:|
| Encoded TIFF / nominal angles | 0.18745 |
| Encoded TIFF / measured global rolls | 0.18506 |
| Encoded TIFF / ideal rear-analyzer transport | 0.18514 |
| ICC gamma decoded / nominal angles | 0.39518 |
| ICC gamma decoded / measured global rolls | 0.38890 |
| ICC gamma decoded / ideal rear-analyzer transport | 0.38914 |

The difference between measured-roll and ideal ray-transport DoLP has median absolute magnitude about 0.00040 within that field. These tests do not remove the broad orientation mismatch. Gamma decoding materially changes magnitude but does not undo the Camera Raw color processing, per-frame white balance, clipping, or look. The decoded values are diagnostic alternatives, not calibrated replacements for the published previews.

Previous museum rendered-TIFF DoLP values, including apparent numerical agreement with a model scenario, remain radiometrically unvalidated. Do not adjust theoretical orientation or ALAN to compensate for this processing issue.

## Local RAW preparation

`src/prepare_spider_linear_raw.py` uses the exact input stems recorded in TIFF XMP: **DSC_5245 through DSC_5256**, corresponding to exposures 1 through 12. It prefers matching NEFs if both NEF and DNG exist. Each original capture timestamp must match the existing TIFF timestamp; exposure duration, ISO, aperture and camera model must remain constant throughout the set.

The converter reads the visible Bayer mosaic, extracts the same physical blue CFA site from each 2×2 Bayer cell, and subtracts the camera-reported black level for that site. It applies no demosaicing, white balance, display gamma, rendered color matrix, Adobe look, per-image brightness scaling or clipping of negative black-subtracted counts. It performs one common linear area resampling to a 2000-pixel long edge, then writes float32 single-channel TIFFs `1.tif`–`12.tif` and `linear_raw_manifest.json`. It is a sensor-channel preparation, not full polarimetric calibration. All original files are read only.

On the computer holding the originals:

```bash
python3 -m venv .venv-raw
source .venv-raw/bin/activate
python -m pip install numpy rawpy exifread tifffile opencv-python-headless
python prepare_spider_linear_raw.py \
  --raw-dir '/path/to/museum/originals' \
  --output museum_linear_raw_01
```

The script is standalone and may be copied from the repository to that computer. Output must be fresh. Original RAWs need not be uploaded. The preparation/mapping and linear Stokes arithmetic have three synthetic unit tests; actual NEF/DNG decoding cannot be validated here because those originals are not available in the current workspace. The new single-channel files require `channel=null` in a new batch registration input. They must not overwrite or be substituted into an existing checksum-locked TIFF batch. New RAW decoding and CFA sampling can change pixel geometry, so a fresh registration is required before north-up/model comparison.

Official API reference for visible RAW mosaic, CFA and black-level access: https://letmaik.github.io/rawpy/api/rawpy.RawPy.html
