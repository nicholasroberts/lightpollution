# Figure 2 spectral and camera-weighting audit — 27 November 2026

Technical working note, not manuscript prose.

## Findings and scope

The camera weighting has a direct physical basis. The selected lamp spectra have a
traceable literature basis as generic source-emission templates. The sky-transfer
step is an explicit approximate scenario; it is not a measured spectral calibration
of either Cottage or museum. Keep that distinction when interpreting the fitted
amounts. No new atmospheric-depolarization parameter is introduced.

Figure 2 should emphasize the spatial DoLP comparison and the required relative
addition at each site. Do not put the LED/HPS percentage in the figure headline or
treat repeated percentages as independent observations. Retain numerical quantities
in machine-readable fit records for reproducibility.

## Camera calibration: checked against the primary report and arrays

ScienceDjinn, Jochen Smolka, Nikon D850 camera calibration report, 24 May 2022,
sections 2 and 6, Figures S6/S7:
https://www.sciencedjinn.com/wp-content/uploads/2022/05/D850-Calibration-Report.pdf

The report measured D850 plus Sigma 8 mm f/3.5 fisheye responses to monochromatic
stimuli. Counts per second per ISO were divided by measured total stimulus PHOTON
radiance. Thus the RGB curves are photon-response curves, not energy-response curves.
For a spectral energy density L_lambda, use integral L_lambda * lambda * R_green /
(h*c) d_lambda. The common h*c and a common channel gain cancel in relative Stokes
ratios. Peak normalization of R_green changes neither fitted physical amounts nor
DoLP when applied consistently to natural and artificial contributions.

Checked data/reference/d850_elf_response.csv against a newly downloaded original:
https://github.com/sciencedjinn/elf/blob/master/calibration/nikon%20d850/rgb_calib.mat
All wavelength samples match every second row of the 0.5 nm source arrays.
Across all RGB channels maximum relative value differences are below 5e-10,
consistent with CSV rounding. No wavelength shift or channel swap found.
ELF classes/Calibrator.m independently documents photon-radiance calibration:
https://github.com/sciencedjinn/elf/blob/master/classes/Calibrator.m

The published camera/lens calibration does not include this project's added
polarizer. Its wavelength-dependent transmission and polarimetric efficiency remain
unmeasured here. A wavelength-independent common transmission cancels; an unknown
spectral transmission does not necessarily cancel between different source spectra.
Do not claim a complete system-specific spectral/polarimetric calibration.

Supporting paper: Nilsson, Smolka & Bok (2022), DOI 10.3389/fevo.2022.951328,
methods and supplementary Figure S1:
https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2022.951328/full

## Lamp templates and transfer: literature basis, conditional sky approximation

src/prepare_spider_spectral_pilot.py selects the author's generic HPS and
LED_2958K-3969K templates from Spectrum Analyser:
https://davinci.fmph.uniba.sk/~kundracik1/SpectrumAnalyser/

Kocifaj, Kundracik & Bily (2020), MNRAS 491, 5586–5594,
DOI 10.1093/mnras/stz3260, section 3.2, describes averaged generic source spectra,
including HPS and three relevant white-LED CCT groups. These are emission
templates, not measured pure skyglow from our sites. The paper discusses relative
spectral power densities and energy radiance units. Applying an energy-to-photon
conversion once in the camera integral is therefore appropriate.
Author-hosted primary paper:
https://davinci.fmph.uniba.sk/~kundracik1/SpectrumAnalyser/stz3260.pdf

The preparation code applies (550/lambda)^4 to the natural source template and
(550/lambda)^exponent to lamp emission, default exponent 2. The latter is a
phenomenological spectral-shape scenario. It is not a complete Rayleigh-plus-aerosol
scattering and extinction calculation, nor a retrieved site atmosphere.
The literature does NOT validate exponent 2 specifically for these observations.

The earlier generated spider_spectral_inputs_01 files were described in the
12 November labbook as a 450–750 nm pilot; later code supports 352–750 nm.
Current Figure 2 takes actual common support within 352–700 nm and camera coverage,
rather than assuming the new generation defaults describe historical files.
Exact EcoVis file provenance must be read, not inferred from the latest source code.
The files are processed local data, not available in this GitHub checkout.

Natural and lamp templates are normalized to equal integrated ENERGY after their
declared transfers over the same actual band. Fitted amounts consequently describe
relative incident sky band radiance. They are not emitted lamp power, lumens,
electrical power, or calibrated absolute sky radiance. Renormalizing after transfer
discards information about absolute propagation efficiency.

Bará, Rigueiro & Lima (2019), JQSRT 239, 106644,
DOI 10.1016/j.jqsrt.2019.106644, sections 2.1/2.3 and Appendix, explicitly combines
source spectra, wavelength-dependent atmospheric propagation and detector sensitivity.
This supports the physical mechanism of detector-dependent responses. It also shows
why no universal lamp-type conversion follows without specifying spectra, geometry,
atmosphere and normalization:
https://arxiv.org/pdf/1909.10909

## Natural template and inference boundary

The current Figure 2 spectral model uses the same transferred natural/Moon template
for Sun and Moon. This preserves the final Figure 1 camera-green mixture under the
separable approximation, but is not a separately calibrated twilight spectrum.
The source Moon spectrum is a source-irradiance input with an approximate sky
transfer, not a measured diffuse spectrum of the capture-time sky.

The fitted unpolarized addition is a conditional explanation of residual DoLP
dilution on the adopted natural baseline. It does not uniquely establish that all
remaining dilution is artificial, or identify the real lamp mixture. Source-specific
scattered light need not be exactly unpolarized; that remains the adopted effective
background approximation. Preserve pmax=1 and no extra depolarization as requested.

For the same shared templates and detector, different sites can require different
relative additions while retaining the same conversion between alternative HPS/LED
amounts. That is physically consistent with the approximation; it is not evidence
that real site skyglow spectra are identical. No model parameters were added merely
to generate differing ratios.

## Actual EcoVis input verification

The new read-only audit checks the actual band, normalization, inherited provenance,
camera support, and hashes/overlaps against BOTH saved v02 fits:

```bash
python src/audit_paper1_fig02_spectra.py \
  --spectra-root data/processed/spider_spectral_inputs_01 \
  --fit-results data/processed/paper1_fig02b_Cottage_02 data/processed/paper1_fig02c_museum_02 \
  --output data/processed/paper1_fig02_spectral_audit_01.json
```

This records rather than invents the original transfer settings and source hashes.
Local validation used synthetic spectra only: actual-band selection, finite overlaps,
normalization, and protection against overwriting a previous audit passed.
No production input spectra or experimental results were regenerated locally.

## Decision

The spectral-overlap mechanism and camera weighting are physically justified.
The exact predicted LED/HPS difference is conditional on approximate sky spectra
and the published camera/lens response, not independently validated for either site.
Continue Figure 2 as a comparison of explicit alternative spectral scenarios,
with the visual spatial match as its principal result. The percentage need not be
quoted. Before making stronger site-specific spectral or biological claims,
evaluate the existing template/transfer and polarizer uncertainties.
