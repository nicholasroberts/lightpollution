# PennCC four-star orientation — 24 November 2026

Investigator selected four stellar points in the locked linear-green reference TIFF and accepted the identified four-star solution. The Moon is not visible. Its pixel position below is a prediction from stars and ephemeris, not an observed Moon centroid. WCpitch's observed candidate [621,638] was separately confirmed.

| Point | Star | Original pixel x | Original pixel y | Geometric elevation | Compass azimuth |
| --- | --- | ---: | ---: | ---: | ---: |
| P1 | Arcturus | 1369.536424 | 634.501538 | 40.908715° | 114.839618° |
| P2 | Capella | 587.864078 | 737.913285 | 35.828734° | 301.122277° |
| P3 | Procyon | 681.578947 | 298.330061 | 26.095602° | 241.670097° |
| P4 | Vega | 1255.263158 | 1120.265975 | 18.760817° | 49.842899° |

Exact clicks and source hashes: config/paper1/PennCC_reference_star_points.json. Reference image is 2000 × 1332 pixels; centre [1001.7451430363744,659.0293225899803], rim radius 655.544999673502 pixels.

## Reproducible protocol

1. Verify selected-point source-lock and reference-TIFF SHA256 checksums. Keep reference centre and radius fixed.
2. Use PennCC NEF DateTimeOriginal 2019:04:20 22:24:20 with explicit +01:00 offset, hence 2019-04-20T21:24:20Z. This accepted date supersedes the older reference table's 22 April for this capture. Site: latitude 52.558887°, longitude −2.145724°, height 0 m.
3. Compute the four identified stars' topocentric alt/az with PyEphem 4.2.1, pressure=0 (no refraction).
4. For original pixel angle beta=atan2(y−cy,x−cx), test A=sign*beta+offset with sign ±1. Estimate offset as the argument of the mean unit phasor exp(i*(A−sign*beta)), equal weight per star. No centre, radius, tilt, source brightness or polarization fitting.
5. Compare full spherical positional residuals using the fixed equisolid image elevation: elevation=90°−2*asin(pixel_radius/(sqrt(2)*rim_radius)). Choose the lower-residual parity for the accepted identities.
6. Best relation: A=−beta+110.620044872646°. This corresponds to production coordinate_matrices handedness='reversed'; verified against the actual north-transform function for all four points. RMS spherical positional residual 1.6763338911610621°. Individual residuals: 2.25476822°, 0.25822155°, 1.30249593°, 2.09600275°. Opposite parity for these same identities gives 53.51047244° RMS. These are registration residuals, not angular uncertainty confidence limits; all four stars are used, with no held-out star.
7. Compute Moon geometric elevation 5.264019536252834°, azimuth 121.55157643192099°, illuminated fraction 97.19844818115234% (waning gibbous). Project using the stellar orientation and fixed equisolid mapping: predicted Moon pixel [1615.1588499211502,540.5543260844313].
8. Investigator noted consistency with the brighter horizon region and accepted use of this prediction for 1D. Horizon brightness is qualitative support, not an additional astrometric datum.
9. Copy source registration trees to a fresh derived root, change celestial-anchor metadata only, and refresh copied lock checksums. Preserve prepared TIFFs, landmark transforms and saved camera Stokes. Store all original hashes and an annotated stellar overlay.
10. Use the predicted Moon as a reference for the existing rotation helper: this encodes the star-derived transform algebraically and does not independently validate the Moon. Suppress the black observed-Moon cross for PennCC in 1D; use only the red predicted marker. Label measured PennCC maps star-anchored.

A prior commentary statement that this stellar solution required the opposite production parity was a label error, caught by direct production-function validation before committing. The implemented convention is reversed, consistent with the existing comparison convention.

## Execution on EcoVis

```bash
git pull --ff-only origin paper-1-polarization-model
python -m pip install -r requirements-stellar-registration.txt

python src/prepare_paper1_stellar_registration.py \
  --registration-root data/processed/spider_registration_green_all_sites_01 \
  --metadata-audit data/processed/nef_metadata_audit.json \
  --output data/processed/paper1_fig01d_registration_01

python src/run_paper1_fig01d.py \
  --registration-root data/processed/paper1_fig01d_registration_01 \
  --raw-manifest-root data/processed/spider_linear_green_all_sites_01 \
  --metadata-audit data/processed/nef_metadata_audit.json \
  --output data/processed/paper1_fig01d_01
```

The earlier 1D attempt stopped in preflight before creating its output directory. Use a new suffix if any output directory already exists.

Preparation exports PennCC_stellar_registration.json/png and derived runs/PennCC and runs/WCpitch locks. The 1D runner consumes per-site declared handedness, calculates signed audited Stokes, transforms I,Q,U, fits one solar fraction per site to spatial AoP changes, and derives both DoLP/AoP from that same mixture. No ALAN or amplitude fit. Linear viridis and romaO remain fixed. Original measured-camera field statistics are distinct from binned/north-resampled diagnostics.

Validation: reproduced the four-star solution and predicted Moon from the uploaded exact clicks; verified every bearing against the actual production coordinate transform. Syntax checks pass. Full real-data preparation, Stokes audit, fit and figure assembly run on EcoVis; no successful field fit is asserted here.
