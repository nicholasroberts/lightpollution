# Figure 1D — PennCC and WCpitch

## Agreed layout

One subfigure D, four rows and four columns:

| Row | Experimental | Lunar only | Solar only | Best combined |
| --- | --- | --- | --- | --- |
| PennCC DoLP | Camera green | Natural model | Natural model | Shared site mixture |
| PennCC AoP | Camera green | Natural model | Natural model | Same mixture as DoLP |
| WCpitch DoLP | Camera green | Natural model | Natural model | Shared site mixture |
| WCpitch AoP | Camera green | Natural model | Natural model | Same mixture as DoLP |

Use linear viridis DoLP 0–1, romaO AoP −90° to +90°, north up/east left, observed and astronomical Moon markers, and the zenith-centred 60° full-field circle.

## Calculation

Reuse the museum fitting method unchanged: fit one solar fraction f per site to spatial AoP changes on training blocks, evaluate lunar-only/solar-only/combined on identical withheld pairs, and report block-deletion stability. Sum (1−f) times lunar I,Q,U and f times solar I,Q,U before deriving either DoLP or AoP. No separate DoLP fraction, ALAN background, dilution, polarization ceiling or sky rotation fit. PennCC and WCpitch can have different f, but each site's two rows must have exactly the same f.

Retain a single training-estimated analyzer-basis offset per site, applied identically to all three theoretical patterns using the existing 1C helper. Report withheld absolute-AoP residuals as well as spatial-change metrics. No independent optical zero calibration is claimed.

## Capture metadata

Investigator confirms PennCC NEF metadata is correct: first retained exposure 20 April 2019 22:24:20 BST (21:24:20 UTC), superseding the 22 April date printed in the older site table for this acquisition. WCpitch first retained exposure is 14 May 2019 22:45:32 BST (21:45:32 UTC). Approximate solar elevations from the metadata audit were −16.9° and −13.5° respectively; the fitting script recomputes exact geometric Astropy solar coordinates from each verified lock.

The runner checks the locked UTC against the first retained NEF timestamp plus explicit offset before any output. It stops on a stale anchor instead of silently changing celestial geometry. Reviewed Moon identity/registration are prerequisites, not inferred from fit quality.

## EcoVis execution

```bash
cd /run/media/bznwr/data/lightpollution
git pull --ff-only origin paper-1-polarization-model

python src/run_paper1_fig01d.py \
  --registration-root data/processed/spider_registration_green_all_sites_01 \
  --raw-manifest-root data/processed/spider_linear_green_all_sites_01 \
  --metadata-audit data/processed/nef_metadata_audit.json \
  --output data/processed/paper1_fig01d_01
```

Existing source locks and files are read only. Uses the existing audited signed nominal Stokes reconstruction, original RAW hash/metadata checks, explicit reversed parity, Moon rotation, and 24° Berry/Wang neutral separation. It does not re-register landmarks.

Outputs per site: audit/, north/, fit/ (solar_lunar_fit.json, arrays, profile, fit diagnostics). Combined figure/: fig01d_two_sites.png/pdf/svg and figure_values.json with fractions, solar/Moon geometry, source hashes, diagnostics and binned central-field summaries. These summaries do not replace original-camera publication statistics.

## Validation and interpretation

Python syntax checked. A synthetic renderer test verified saved-mixture consistency, reuse of the 1C analyzer-basis alignment, derived combined DoLP and rejection of tampered mixture arrays. The test used stub geometry components; it does not validate the physical sky model or the actual site measurements. EcoVis execution and review remain required. Report the actual fit outcomes, including boundary fractions or no improvement, rather than presupposing a better combined fit.

## Accepted PennCC stellar anchor — supersedes observed-Moon prerequisite for PennCC

PennCC now uses the investigator-accepted four-star orientation, with Moon pixel predicted from that orientation and ephemeris, not observed. WCpitch Moon [621,638] is confirmed. Run the new derived-lock preparation before the 1D runner; use its output as --registration-root. Full protocol and commands: [PennCC stellar orientation](PennCC_stellar_orientation.md). Both sites use the production reversed parity; this is directly tested from stars for PennCC. Suppress its observed-Moon cross and retain the predicted marker. Each site still has one shared solar fraction across DoLP/AoP.
