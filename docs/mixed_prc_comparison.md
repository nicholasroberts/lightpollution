# Camera-band PRC comparison after analyzer-sign checking

`src/compare_spider_mixed_prc.py` compares measured green-camera Stokes with
the saved lunar-only, solar-only and mixed fits. It does not refit source
weights, background, geometry or analyzer zero. Lock and Stokes checksums and
all saved model arrays are checked before calculation. The explicit analyzer
convention recorded by `correct_spider_analyzer_sign.py` is required.

The sampled field is the intersection of the central 60-degree full field,
the saved fitting support, and physically valid Stokes for every model and the
measurement. Every comparison uses identical included pixels. Coverage is
reported; excluded pixels are not filled. Equal binned image pixels are used
under the existing provisional equisolid projection. Bin size must match the
saved fit (default 8).

Sum I, Q and U separately. With parallel/perpendicular absorbances 0.4/0.04,
let b = 1 - (10^-0.4 + 10^-0.04)/2 and a = (10^-0.04 - 10^-0.4)/2.
Then R1 = b sum(I) + a [sum(Q) cos(2 theta) + sum(U) sin(2 theta)] and R2
has the opposite modulation. Signed PRC is ln(R1/R2).

Writing k = a hypot(sum(Q),sum(U))/(b sum(I)), the maximum absolute PRC is
2 atanh(k). The median absolute PRC over a continuous uniform distribution
of receptor orientations is 2 atanh(k/sqrt(2)). Both are unchanged by a
constant rotation of the Stokes zero. Absolute compass-heading agreement is
not claimed. Native-axis curves at 1-degree steps are saved separately.
Displayed peak-relative curves have each pattern's own mathematical peak as
zero; this is an amplitude comparison, not a fitted angular correction.

This is the existing camera-plane integrated receptor proxy, not yet a
calibrated retinal/optical transport or pigment-spectral prediction. It uses
actual measured intensity weighting for the measurement, and the saved
DoLP-fit proxy intensity weighting for each model. The uniform fitted
background retains its phenomenological interpretation. Comparison over
the fitted support is a derived check rather than independent validation.
The additional withheld-block calculation samples a fragmented field and
must not be interpreted as a second full-eye result or independent night.

## Run on EcoVis

```bash
git pull origin publication-figures

for site in WClawn oxfordbussinesspark; do
  python src/compare_spider_mixed_prc.py \
    --reference-lock data/processed/spider_sign_checked_01/registration/locked_site_references.json \
    --fit "data/processed/spider_sign_checked_01/$site/fit" \
    --output "data/processed/spider_prc_checked_01/$site"
done

python -m zipfile -c data/processed/spider_prc_checked_01_results.zip \
  data/processed/spider_prc_checked_01
```

Each site produces `prc_comparison.json`, `prc_comparison.png`,
`prc_orientation_curves.csv` and `prc_support.npz`. Use a fresh output folder.
