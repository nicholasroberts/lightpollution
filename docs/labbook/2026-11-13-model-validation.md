# Fixed-geometry and natural-envelope audit

This adds steps 1 and 2 of the model validation, without modifying the saved
measurements, registration, historical PRC or spectral scenario inputs.

## Run on EcoVis

From the lightpollution repository on publication-figures:

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/validate_spider_sky_model.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
  --neutral-separation-deg 10 \
  --output data/processed/spider_model_validation_01
```

The 10 degree neutral-point separation is an explicit scenario assumption,
not a fitted value. Use the same separation as the comparison being audited;
if changing it, use a separate output directory. Existing outputs are protected.

Look in **data/processed/spider_model_validation_01**:

- SITE_validation.png: measured DoLP, undiluted natural envelope, raw AoP
  residual and mean DoLP against angular separation from the Moon.
- validation.json: raw and offset-removed axial AoP residuals/coherence;
  measured-to-envelope ratio quantiles, means, medians, errors and fraction
  exceeding the envelope for atmospheric factors 1 and 0.8; zenith 30 degree
  radius and extended 60 degree radius regions are separate. The latter is a
  diagnostic region, NOT the animal's 60 degree full field.
- site_inventory.csv: all 18 sites, candidate low-ALAN status, pending anchors
  or missing inputs. Only sites with anchors are evaluated. Initially these
  are Cottage and Oxford; other sites need independent Moon review first.

## Interpretation and safeguards

Moon elevation is fixed by the saved ephemeris; image bearing is fixed by the
reviewed candidate. No source search or amplitude optimization occurs. The
candidate's measured equisolid elevation minus ephemeris elevation is reported
as a camera-centre/radius/tilt diagnostic, not corrected away.

The existing global image-plane Stokes basis is retained. Absolute AoP zero,
handedness and ray-dependent instrument response are not established. Raw AoP
residuals are therefore convention diagnostics. A per-site constant offset is
also reported to describe shape agreement; this is explicitly NOT independent
validation of absolute AoP. No offset-removed map is presented as a prediction.
A future calibration must fix the convention independently and then assess
held-out sites.

The low-ALAN candidate screen uses the historical reference mean VIIRS radiance
<=1, independent of measured DoLP. It does not prove an ALAN-free sky, nor is
VIIRS a ground skyglow calibration. Cottage, EastHorrington, Gavarnie and
slabhouse qualify under this declared screen. Missing Moon anchors stay pending.
Neither the Cottage target DoLP=0.45 nor historical PRC enters this audit.
Ratios to the natural envelope cannot distinguish atmospheric depolarization,
ALAN, spectral response, geometric calibration or model errors by themselves.
The Berry-normalized magnitude as a natural DoLP envelope remains an assumption.

Terrain is not automatically identified from DoLP. Default regions are angular
cuts and remain provisional: trees can enter even a zenith field. To supply
reviewed masks, use --sky-masks PATH, containing SITE.npy boolean arrays in the
saved Stokes array coordinates (True=sky). Masks must match shape exactly and
are hashed. Do not use a high-DoLP threshold to define sky: that biases the
natural-envelope check. Both measurements and models must be finite and physical
for comparisons; invalid/noncomparable pixels are counted, never clipped.
AoP additionally requires both DoLP values >=0.1; this threshold is recorded.
Equisolid pixels have equal solid angle (the code assumes that lens mapping).
The historical camera channel is not an independently measured green510 spider
response, so this is not a species-specific spectral validation.

## Verification

Four numerical tests cover axial wraparound, retaining an imposed AoP offset,
recovery of a known descriptive envelope ratio, and empty eligible samples.
A synthetic saved-registration run checks inventory, JSON and figure generation.
Real sky results require running the command on EcoVis; none are claimed here.
