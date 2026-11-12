# Wang 2016 dynamic spider model

Source: Wang, Gao, Fan & Roberts (2016), Journal of Optics 18 065601, DOI 10.1088/2040-8978/18/6/065601. Investigator-supplied PDF SHA256 70aebbc22e9855fc6c127933b7456148b782fbc56bec9bec1b88dbf60dd15429. Equations inspected in rendered PDF, not text extraction alone.

## Model and implementation

`src/wang2016_sky.py` implements printed complex-field eq1/3, neutral points eq4/5 with delta=4atan(L), turbidity eq6/7 exp(-T/4+.12), horizon eq8 cos(theta)^(1/10), Perez eq12 with caller-supplied coefficients, spectral basis eq13-15 and S eq16, and full combined eq17 when S is explicitly supplied. Dimensionless stereographic coordinates z=tan(theta/2)exp(i phi) are distinct from the fisheye image projection and from the display scale in eq18. AoP=arg(w)/2 in the paper's stated convention. The arbitrary-azimuth neutral-point equations are transcribed as printed, not silently replaced by a different geometric model; their coordinate/basis correspondence with camera Stokes remains an empirical verification task. Root-at-zero analytic limits are rejected explicitly. Nonphysical Stokes outputs are rejected, never clipped to create agreement.

The full eq17 spectral branch needs verified daylight basis tables and the supplementary coefficient derivations S1-S5. These are not in the supplied main PDF, and the publisher supplement could not be retrieved in this session. No coefficients are invented. The paper's spectral sunlight construction also needs a justified lunar adaptation before making HPS/LED conclusions. Relative monochromatic pilot is clearly separate from full spectral reconstruction.

## Additive ALAN

I_total=I_natural+I_ALAN, Q_total=Q_natural, U_total=U_natural. Hence p_total=p_natural I_natural/(I_natural+I_ALAN), AoP unchanged where defined. ALAN may be an array, allowing spatial structure. Uniform relative dilution is equivalent at the DoP-amplitude level to deltaT=4 ln(1+I_ALAN/I_natural) ONLY when this ratio is constant across directions. It does not make aerosol extinction/phase functions, Perez chromaticity changes or other atmospheric effects identical to ALAN.

`--replace-turbidity-with-alan` sets the natural amplitude multiplier to 1 and adds explicit ALAN, rather than counting both sources of depolarization. This does not imply T=0 in the paper's eq7, which gives exp(.12)>1. The pilot implements eq6 times eq8 as a declared reduction, not full eq17.

## First runnable check

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin publication-figures
source .venv/bin/activate
python src/run_wang2016_spider_pilot.py \
  --demo \
  --neutral-separation-deg 24 \
  --replace-turbidity-with-alan \
  --output data/processed/wang_spider_angular_pilot_01
```

24deg separation is an EXPLICIT ILLUSTRATIVE input, not a recovered constant or site fit. Demo sweeps source elevation 5..85deg at fixed azimuth 30deg; it is not an actual Moon track. Natural relative intensity is uniformly 1 at every frame; ALAN additions 0,.5,1,3. Animal field is zenith-centred, 30deg radius, equal equisolid image-pixel contributions. Eye absorbances .4/.04; signed log ratio versus body heading is exported, with maximum absolute contrast as an envelope. No fitted VIIRS curves, historical observations, image alignments or measured datasets changed.

Outputs: wang_spider_pilot.png, heading_responses.csv, provenance.json. Upload PNG and provenance first. Small pilot completed successfully here both with retained T3 amplitude and with amplitude replaced by explicit ALAN. Eight targeted tests pass. These are analytic/synthetic checks, not empirical validation.

## Actual lunar track

Install `python -m pip install -r requirements-wang.txt`, then replace --demo with --start-utc and --end-utc (explicit Z or +00:00), --latitude and --longitude. The code uses Astropy get_body('moon') with builtin ephemeris, explicit site/UTC and pressure=0; network IERS updating disabled and relevant precision warnings retained. This geometric path was checked against official API documentation but was NOT executed in the current workspace (Astropy unavailable). No claim to exact high-precision ephemeris or atmospheric refraction. Dates/times must be independently verified before comparing to the original photograph. Earlier TIFF inspection records can contain export/clock anomalies; do not treat Photoshop DateTime as automatically verified exposure UTC. Maximum interval two days, positive step required; source below horizon skipped.

Moon phase, source irradiance, lunar spectrum, changing atmospheric extinction, twilight/solar contribution and natural background are not simulated in this reduced pilot. It supplies geometric motion under a controlled intensity assumption. Solar model extension to lunar positions is provisional; full wavelength-dependent lunar validation remains outstanding.

## Next empirical step

Recover supplement/original author implementation and verify eq17 spectral branch; confirm camera time/zone and field-date metadata; establish camera azimuth from Moon position; compare predicted DoP/AoP at capture time with sky-only image data. Trees mask the celestial-pattern fit but must remain obstructions in the subsequent animal exposure calculation. Do not fit the model to historical PRC alone or infer an HPS/LED ratio from the illustrative intensity sweep.
