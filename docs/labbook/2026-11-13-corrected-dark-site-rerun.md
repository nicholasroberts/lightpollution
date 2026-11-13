# Correct visible Moon bearings and freeze the camera AoP convention

The uploaded real-data audit favours conjugation: Cottage MAE 3.545deg versus
15.437deg direct; held-out Oxford 5.399deg versus 30.692deg direct. All Oxford
offsets came from Cottage unchanged. The selected camera comparison is
conj(W)*exp(2i*offset), offset=-80.37606797554504deg in this run. Conjugation
reverses U and preserves DoLP. It does NOT justify conjugating the physical
animal-response model in a different coordinate basis. No production animal
model basis is changed here. The residual Oxford mean axial shift is -4.469deg.

The automatic source candidates for the new sites hit edges. Explicit visual
bearing estimates, read from the investigator-uploaded intensity overlay, are
now in config/spider_dark_site_moon_bearings_20261113.json:

| Site | Rejected automatic x,y | Approximate visible source x,y |
|---|---|---|
| Skillinge | 637,193 | 650,1135 |
| slabhouse | 1517,440 | 1570,600 |

These are approximate raw-pixel positions, with roughly +/-25px visual-reading
uncertainty, not measured centroids or statistical intervals. They set bearing
only. Moon elevations remain exactly the saved ephemeris values 12.5925536185deg
and 11.1361209834deg; timestamps and registration hashes are preserved. No
polarization map or PRC fit chooses the bearing. The corrected review figure
shows visible-source estimates, the ephemeris-projected model sources and old
rejected candidates. Inspect it before interpreting the new sites.

## Run on EcoVis

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/run_spider_dark_site_validation.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --moon-anchors \
    data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
    data/processed/moon_skillinge_slabhouse_review_01/moon_ephemeris_anchors.json \
  --fixed-offset-audit data/processed/spider_aop_dark_audit_01/convention_and_atmosphere_audit.json \
  --output data/processed/spider_dark_site_validation_01
```

This consumes the two anchor files and previous AoP audit already produced.
The output directory must be fresh. No raw TIFF timestamp read, registration
search, Moon ephemeris recalculation or new AoP offset fit occurs. The neutral
point separation defaults to the previous 10deg scenario. Prior audit files
that record it must agree; older files did not include that field.

Look in **data/processed/spider_dark_site_validation_01**:

- corrected_moon_review.png: check the two visual bearing estimates.
- three_site_dop_comparison.png: original measured DoLP at Cottage, Skillinge
  and Slabhouse against prescribed atmospheric factors 1.0, 0.8, 0.6 and 0.4.
- summary.json: selected conjugate convention residual metrics for all four
  sites, unchanged common offset, radiances and input hashes.
- audit/: per-site six-convention residual panels, DoLP profiles, CSV metrics,
  inventory and detailed JSON, with all prior offsets frozen.

All four sites must produce diagnostics or the runner reports an incomplete
run. A partial output is retained for diagnosis, not called complete.

Atmospheric factors are scenarios, not fitted values. Without independently
justified capture-camera-band ALAN/natural ratios, plots show natural-only
predictions and explicitly leave ALAN unspecified. Real site radiances remain
Cottage 0.36380952381, Skillinge 5.555, Slabhouse 0.498235294118, Oxford
41.8296153846 nW cm^-2 sr^-1. These historical upward radiances are not silently
converted to ground ALAN/natural ratios. The Cottage 0.45 scenario target is
not used. Optional --capture-alan-ratios and --sky-masks pass through to the audit.

Terrain masks remain absent by default; nominal equisolid geometry, exposure
calibration and camera spectral response limitations remain. The workflow is
prepared and tested locally with synthetic Stokes; actual four-site numerical
results require this EcoVis run. No actual new-site validation is claimed yet.

## Verification

Seven unit tests cover handedness/held-out offset transfer, spin-2 terms,
atmosphere-versus-ALAN separation, invalid ratios, preserving anchor elevation,
UTC and report hash when replacing bearing, and conflicting/missing anchors.
A complete four-site synthetic runner exercises the frozen-offset subprocess,
corrected source review, combined DoLP plot, inventories and summary JSON.
