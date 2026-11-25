# Figure 2 — spectral ALAN addition to the accepted natural patterns

Updated25November2026; replaces the earlier PRC-matched Figure2 sketch.
Approved Figure1 remains unchanged.

| Panel | Content |
| --- | --- |
| 2A | Moon/natural, HPS and LED input spectral shapes; equal integrated radiant energy on a common wavelength band |
| 2B | Cottage measured DoLP versus independently DoLP-fitted HPS-only and LED-only predictions using final Figure1 natural fit; report required amounts and LED/HPS percentage difference |
| 2C | Museum, same comparison and fit definition |
| 2D | 510nm pigment photon sensitivity and published ELF D850 green photon response |
| 2E | Predicted PRC for the best DoLP-fitted HPS/LED scenarios at both sites, with experimental proxy where specified; no second fit to PRC |

## Fixed scientific settings

Use unity-amplitude natural Berry/Wang model (no imposed pmax below1),
existing horizon correction N=10, NO added atmospheric depolarization.
Retain source ephemerides, accepted geometry, analyzer conventions and each
site's exact final Figure1 solar fraction and shared neutral separation:
Cottage1.1626536%solar,69.2184859deg; museum74.3078968%solar,70.7216738deg.
Added ALAN is unpolarized and initially uniform, contributing to I only.
HPS and LED are separate alternative fits, not an inferred mixture of site lamps.
No atmospheric parameter is fit alongside the lamp amount.

2B/2C objective is measured camera-green DoLP, explicitly superseding the older
PRC-only Figure2 fit proposal. Specify common60degree field support and exact
DoLP objective before implementing those panels; assess maps/residuals as well.
Use camera response in the forward comparison; then calculate PRC for the SAME
physical fitted lamp amounts with the animal-response approximation.
Do not refit lamp amounts when plotting PRC. The camera-to-animal approximation
and inherited common Sun/Moon spectral template remain declared.

## 2A implementation

src/plot_paper1_fig02a.py reads exactly the HPS_visible_input.csv,
LED_mid_visible_input.csv and spectral_provenance.json in the selected spectral
input directory. Default data/processed/spider_spectral_inputs_01 avoids the
known conflicting alternate pilot copies. Require valid monotonically
increasing wavelength axes and finite nonnegative spectra; verify the two
natural templates agree after normalization. Use existing transferred natural
column as Moon template and alan columns as lamp inputs. Retain native spectral
sample wavelengths plus1nm shared grid, with no extrapolation. Common band
is restricted to352–700nm and actual available support. Normalize each energy
spectrum to integral1 over that band. No camera/animal/human weighting in2A;
no fitting and no extra atmosphere applied. Preserve inherited spectral provenance
inside output JSON, with source hashes. No new site-specific lunar spectrum implied.

Render oneA label, Moon blue, HPS orange, LED green; linear spectral-energy
axis. Output PNG/PDF/SVG, normalized spectral CSV, full provenance JSON and a
verified versioned figures-folder copy. No production spectra are fabricated
or downloaded. Exact lamp/source provenance remains that of the chosen inputs.

```bash
python src/plot_paper1_fig02a.py \
  --spectra-root data/processed/spider_spectral_inputs_01 \
  --output data/processed/paper1_fig02a_spectra_01 \
  --figure-copy figures/paper1/fig02a/v01
```

2A input/normalization tests use synthetic spectra only. Actual input checks,
render and snapshot execute on EcoVis. 2B–2E are planned, not implemented by2A.

Theory: [shared framework](../../theory/theory.md).
