# Figure 2 — spectral ALAN addition to the accepted natural patterns

## Jones spectra adopted for production — 27 November 2026

The investigator accepted the Jones comparison and requested production figures.
Run src/run_paper1_fig02_jones.py to make new Figure 2A–C v04 snapshots using
the saved official Jones/ESO FITS, without repeating online requests.

Jones replaces the LUNAR spectral template only. Retain the previous natural
source-times-lambda^-4 template as the existing approximate SOLAR proxy; there is
no independently modelled twilight spectrum. Do not label the lunar model a
solar prediction. The Figure 2A plot explicitly distinguishes both lunar curves,
the retained solar proxy, HPS and LED.

Keep the final Figure 1 camera-green mixture f exactly fixed. With photon-energy
overlaps K_M and K_S, calculate the band-energy weights by normalizing

    w_M = (1-f)/K_M
    w_S = f/K_S

so the integrated camera-green I/Q/U recover the locked pattern, up to one common
intensity scale. These are deterministic conversions, not fitted parameters.
Save both camera-green f and band-energy weights; they are different quantities.
For Jones runs the saved natural_I/Q/U are band-energy component sums, and the
saved camera_natural_I/Q/U are response-integrated sums. Separate component arrays
and spectra are retained for later receptor calculations.

Use the site-specific Jones zenith spectral shape uniformly over the field as
the explicitly simple spectral approximation. Retain the fixed ESO reference
atmosphere, 400–700 nm band, D850 response, generic lamp templates and ALAN
transfer exponent. No extra wavelength^-4 or extinction is applied to Jones output.
Fit HPS and LED independently to spatial camera-green DoLP exactly as before.
The spectral update changes inferred band-energy amounts, not the approved
camera-green spatial pattern or fit objective. Do not claim improved spatial
residuals simply from changing this scalar spectral conversion.

The runner checks source Figure 1 fit hashes and raw Jones FITS hashes before
using the predictions. It preserves all prior artwork, makes v04 for A/B/C,
and omits percentage claims from figure headlines and its combined summary.
Quantitative per-fit records remain reproducible.

```bash
git pull --ff-only origin paper-1-polarization-model
python src/run_paper1_fig02_jones.py \
  --jones-root data/processed/paper1_jones_moon_spectra_01 \
  --spectra-root data/processed/paper1_spectral_inputs_02 \
  --output data/processed/paper1_fig02_jones_01 \
  --figure-version v04
```

Dependencies are those already installed for the Jones comparison.
Source arrays and aggregate fig02_jones_summary.json stay in the processed output;
artwork and per-panel provenance are copied to figures/paper1/fig02a/v04,
fig02b/v04 and fig02c/v04.


## Approved visible band — 27 November 2026

Use **400–700 nm** consistently for Figure 2A/B/C and the subsequent
camera/animal comparison. This supersedes the historical commands/band choices
below. The 450 nm lower edge in existing v01/v02 fits came from the older
NIST2013 lunar-table pilot, not from the D850 response.

The shared spectral loader now requires complete input coverage of 400–700 nm.
It fails on 450-start input files instead of narrowing the integration band.
Each spectrum is normalized to unit radiant energy over exactly 400–700 nm;
energy-to-photon weighting remains wavelength times photon response.
Source tables can cover a wider interval, but only the approved band is used.

Generate fresh wider-coverage inputs using the existing audited QASUME/air-LUSI
lunar preparation and generic lamp source templates. This replaces the coarse
historical NIST lunar shape as well as extending the band; changed fitted
amounts should not be attributed solely to inclusion of 400–450 nm.
Retain natural transfer exponent 4 and ALAN exponent 2. No missing spectral
range is extrapolated from the old 450 nm file. No UV receptor is added.
The current composite-source overlap check remains a source-shape check,
not site-specific spectral validation.

All Figure 1 geometry, proportions and separation remain fixed. Preserve old
artwork and fits. New Figure 2A uses v02; new B/C use v03. The shared figure
headline no longer quotes the LED/HPS percentage.

```bash
git pull --ff-only origin paper-1-polarization-model
python -m pip install -r requirements-spectral.txt

python src/prepare_spider_spectral_pilot.py \
  --alan-transfer-exponent 2 \
  --output data/processed/paper1_spectral_inputs_02

python src/plot_paper1_fig02a.py \
  --spectra-root data/processed/paper1_spectral_inputs_02 \
  --output data/processed/paper1_fig02a_spectra_02 \
  --figure-copy figures/paper1/fig02a/v02

python src/fit_paper1_fig02b.py \
  --spectra-root data/processed/paper1_spectral_inputs_02 \
  --output data/processed/paper1_fig02b_Cottage_03 \
  --figure-copy figures/paper1/fig02b/v03

python src/fit_paper1_fig02c.py \
  --spectra-root data/processed/paper1_spectral_inputs_02 \
  --output data/processed/paper1_fig02c_museum_03 \
  --figure-copy figures/paper1/fig02c/v03
```

Preparation may retrieve the author lamp archive and lunar source datasets;
requirements-spectral.txt includes the existing spreadsheet/NetCDF readers.
No RAW image processing or Figure 1 refitting is needed.


Updated 27 November 2026; replaces the earlier PRC-matched Figure2 sketch.
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
PRC-only Figure2 fit proposal. 2B implements common physical60degree support on the saved binned Stokes grid and independently minimizes equal-pixel spatial DoLP squared residuals for each lamp; report RMS, mean residual and optimizer diagnostics. Apply this same definition to2C.
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


## 2B implemented — Cottage

src/fit_paper1_fig02b.py retains final1C natural component arrays and independently fits uniform unpolarized HPS/LED additions to the full common-field spatial camera-green DoLP pattern using ELF D850 green photon sensitivity. See25November labbook for precise normalization/support. Four columns measured/natural/HPS/LED with linearviridis0–1. Each lamp has its own nonnegative optimizer call. The observed mean is reported but is not constrained to match. Fixed scalar spectral overlaps still imply equivalent camera predictions and a site-independent LED/HPS amount ratio under this model. Report spatial residuals and store full fit/spectral/input provenance. Fitted quantities are relative to total natural sky band radiance, not lunar-only or electrical units. Preserve same physical lamp amounts for subsequent animalPRC in2E. No extra atmosphere or natural geometry refit.

Command: `python src/fit_paper1_fig02b.py --spectra-root data/processed/spider_spectral_inputs_01 --output data/processed/paper1_fig02b_Cottage_02 --figure-copy figures/paper1/fig02b/v02`.


## 2C implemented — museum

src/fit_paper1_fig02c.py uses the identical shared2B calculation with museum final1C data, panelC and museum output names. All objective, spectra, D850 response, normalization and saved-array conventions are shared; no new atmospheric or natural-pattern fit. Output figures/paper1/fig02c/v02/fig02c_museum.png with report lamp_DoLP_fit.json. Figures2D/2E remain next work.

The v01 outputs used field-mean matching. Preserve them as historical outputs; use fresh _02 data directories and v02 figure snapshots for the spatial fits. Figure 1 natural fits and approved artwork are unchanged.
