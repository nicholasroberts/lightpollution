# PRC intensity-weighting diagnostic

`audit_spider_prc_weighting.py` reuses the signed-lock PRC comparison to verify
checksums, saved models and identical central-field support. The `baseline/`
subfolder reproduces the preceding PRC results. No model parameters or
analyzer zero are adjusted.

The diagnostic combines each intensity-weight map with each complex
fractional polarization field (Q+iU)/I. Matrix rows provide intensity weights;
columns provide fractional polarization. Diagonal entries reproduce original
PRC; the measured-intensity row compares every pattern under common observed
weighting, and the uniform row removes spatial intensity variation.
These swaps are counterfactual diagnostics, not corrected physical models.

Integrated DoLP factors exactly into:

1. Intensity-weighted mean local DoLP: sum(I*DoLP)/sum(I).
2. AoP cancellation factor: abs(sum(I*DoLP*exp(2i*AoP)))/sum(I*DoLP).

The cancellation factor is one for aligned polarization and zero for complete
cancellation. It is undefined (JSON null) if polarized signal is identically
zero. Both factors and the resulting PRC amplitude are invariant to a constant
analyzer-zero rotation. Contributions interact; the factorization does not
assign unique additive causal error shares.

Intensity maps are normalized by their own mean on the common support, with
normalized residual statistics. This checks spatial distributions only, not
absolute radiometry. Full-field and fragmented withheld-block calculations
remain separate. The inherited camera-band and camera-plane limitations apply.

```bash
cd /run/media/bznwr/data/lightpollution
git pull origin publication-figures

for site in WClawn oxfordbussinesspark; do
  python src/audit_spider_prc_weighting.py \
    --reference-lock data/processed/spider_sign_checked_01/registration/locked_site_references.json \
    --fit "data/processed/spider_sign_checked_01/$site/fit" \
    --output "data/processed/spider_prc_weighting_01/$site"
done

python -m zipfile -c data/processed/spider_prc_weighting_01_results.zip \
  data/processed/spider_prc_weighting_01
```

Use a fresh output directory. Outputs include `prc_weighting_audit.png` and
`prc_weighting_audit.json`, plus the unchanged-baseline comparison products.
