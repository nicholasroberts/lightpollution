# AoP convention and dark-site DoLP checks

The investigator requested Cottage, Skillinge and Slabhouse as checks on the
natural DoLP model, retaining actual site radiances and allowing atmospheric
depolarization. These are NOT identically dark sites:

| Site | Reference mean radiance (nW cm^-2 sr^-1) | Ratio to Cottage |
|---|---:|---:|
| Cottage | 0.36380952381 | 1.000 |
| Skillinge | 5.555 | 15.269 |
| slabhouse | 0.498235294118 | 1.369 |
| oxfordbussinesspark | 41.8296153846 | 114.977 |

Exact means, SDs and ratios are read from prc_spider_18_sites_verified.csv in
every run and copied to the output inventory. Values above are rounded context.
The input is historical multiyear upward radiance, not measured capture-time
ground skyglow. It cannot directly supply an ALAN/natural intensity ratio.

## First run: Cottage calibration, Oxford check

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/audit_spider_aop_conventions.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
  --neutral-separation-deg 10 \
  --output data/processed/spider_aop_dark_audit_01
```

Look in that single output directory. SITE_aop_conventions.png compares six
explicit spin-2 hypotheses: direct/conjugated phase, each with zero or +/-2phi
basis term. Each hypothesis gets ONE angle offset from Cottage. Oxford and any
additional sites use that offset unchanged. aop_metrics.csv and
convention_and_atmosphere_audit.json report all alternatives, without selecting
or changing the production convention. A good result must transfer to held-out
sites; a better Cottage fit alone does not establish the convention. These
alternatives diagnose handedness/global-versus-radial basis errors, not a full
ray-dependent optical calibration. They do not move source coordinates or fit
spatial parity. Camera exposure/calibration errors remain possible.

Cottage_dop_atmosphere.png shows original measured camera DoLP against natural
eta=1, 0.8 and 0.6 profiles. No Cottage target=0.45, historical PRC fit, or fitted
DoLP scale enters. The white circle in earlier maps means a 30deg radius field;
this audit reports zenith-field DoLP ratios separately from 60deg-radius AoP
checks. Both exclude 15deg around the Moon. AoP requires measured/model DoLP>=0.1.
Optional --sky-masks DIR reads SITE.npy boolean sky masks in saved Stokes array
coordinates. Without these, every result remains provisional for terrain.

## Add the other two dark-site checks

They need their own ephemeris and independently reviewed image bearing. Prepare
review panels from their original TIFF capture metadata:

```bash
python src/anchor_moon_from_local_time.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --sites Skillinge slabhouse \
  --output data/processed/moon_skillinge_slabhouse_review_01
```

The new --sites option preserves the old Cottage/Oxford defaults. Sweden uses
Europe/Stockholm, Slabhouse Europe/London, with date-specific daylight saving.
Inspect timestamps, reference-date flags, cyan elevation rings and candidate
locations in moon_ephemeris_anchors.png. A TIFF DateTime fallback may be an
export date. Never approve a candidate simply because it fits a polarization
pattern. --local-time-overrides accepts a JSON object of independently confirmed
SITE: "YYYY:MM:DD HH:MM:SS" timestamps when necessary. Raw images are unchanged.

Only after independently checking the new candidates, run:

```bash
python src/audit_spider_aop_conventions.py \
  --batch data/processed/spider_prc_all_sites_01 \
  --moon-anchors \
    data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json \
    data/processed/moon_skillinge_slabhouse_review_01/moon_ephemeris_anchors.json \
  --reviewed-sites Cottage oxfordbussinesspark Skillinge slabhouse \
  --neutral-separation-deg 10 \
  --output data/processed/spider_aop_dark_audit_02
```

This adds all three dark-site DoLP profiles and held-out AoP checks.
site_radiance_inventory.csv makes pending sites explicit in the first run.
The 10deg neutral-point separation remains a prescribed model assumption.

## Atmosphere versus ALAN

For a declared uniform capture-band addition a=I_ALAN/I_natural,
p_measured / p_natural = eta/(1+a). DoLP alone cannot identify eta and a.
Without independent a, the JSON therefore reports *apparent attenuation* only.
It does not call that ratio atmospheric depolarization. The plotted 0.6/0.8
factors are sensitivity scenarios, not fitted atmospheric values.

Optional --capture-alan-ratios PATH accepts JSON {SITE: ratio}, with nonnegative
finite camera-band capture ratios justified independently. With these supplied,
the script also reports conditional eta quantiles and the fraction exceeding 1.
It does not clip those estimates or optimize them. Do not substitute raw VIIRS
radiance for this dimensionless ratio or silently use spider spectral ratios for
the historical camera channel. Full HPS/LED, moon-brightness and transfer
assumptions remain in the separate spectral forward workflow; this audit does
not overwrite their configuration.

## Checks completed here

Four tests check conjugation with a common calibration offset transferred to a
held-out synthetic sky, spin-2 basis changes, separation of atmospheric and ALAN
terms, and invalid ratios. End-to-end two-site synthetic figures/JSON recover
the imposed conjugated convention on the held-out site. Existing Moon spectrum
checks also pass. These synthetic results are software checks, not evidence
about the real skies. Real comparisons require the saved EcoVis measurements.
