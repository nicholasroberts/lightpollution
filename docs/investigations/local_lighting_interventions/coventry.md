# Stage A case dossier — Coventry

**Original NUTS2016 container:** UKG33 — Coventry

**Mechanism:** Operational switch-off before principal citywide retrofit

**Intervention period:** May–July 2024 staggered switch-offs

**Status:** Primary documentary evidence found, but the precise treated lamp-by-lamp historical polygon and suitable internal controls are **not yet certified**. No VIIRS or ISS local result has been measured.

## Scientific test

Test a 2024 DNB radiance reduction from hours-of-operation policy largely preceding the citywide 2025 LED contract; compare switched-off eligible streets with independently verified exemptions. No inference that every 2024 lamp was non-LED.

## Independently grounded treatment footprint and comparison

**New tangible source:** council publishes 51-page *Part Night Switch Off List of Streets (and Ward)* and downloadable XLSX. Use *street × ward* as a provisional affected unit. Some long streets may be partially treated, and a listed street does not identify precise lamp columns or commissioning date. Compare with **confirmed excluded** city-centre streets, hazard locations or high-risk road lamps **only after** documented CMS exemptions; 'absent from list' does **not** automatically mean control. Seek per-column asset GIS and May–July 2024 programming/activation log.

## Dated documentary evidence and caveats

20 Feb 2024 budget vote; CMS programming began May 2024 and all eligible lamps were programmed by July. Approx **60%** included after exemptions. Weekdays generally 00:00–05:30; weekends 01:00–05:30. Council FOI response dated 8 Oct 2024 gives **31,400 units / 32,300 lanterns** in PFI stock, reflecting multi-head columns; avoid confusing lamps and columns. Main LED replacement rollout started **Nov 2025** and restoration of 2024 overnight lighting was in **Sep 2025**, both outside locked 2014–24 study.

## Records required to lock Stage A footprint

Download **XLSX** and PDF of switched-off streets; retain published version/date and map streets to OS/OpenStreetMap centre-lines with spatial review, not automatic pixel classification. Request 2024 CMS log with lamp ID, exact blackout dates/times, latitude/longitude, type (HPS/LED) as of May 2024, exceptions, and pre/post electricity. Need 2023–24 monthly VIIRS; first clean post-programming months Aug–Dec 2024; 2024 annual product is mixed-exposure.

## Original primary reports and corroborating documents

1. [Council — Part Night Switch Off streets/wards (PDF 51pp + XLSX)](https://www.coventry.gov.uk/downloads/download/7837/part-night-switch-off-alphabetical-list-of-streets-and-ward) — source street register; date of publication not explicit; available 2026.
2. [Council — source street PDF](https://www.coventry.gov.uk/downloads/file/43090/part-night-switch-off-list-of-streets-and-ward-pdf) — PDF 51pp; undated.
3. [Council — source street XLSX](https://www.coventry.gov.uk/downloads/file/43091/part-night-switch-off-list-of-streets-and-ward-excel) — original XLSX; undated.
4. [Council — Jan 2025 scrutiny report/minutes](https://edemocracy.coventry.gov.uk/mgAi.aspx?ID=43074) — primary record for 2024 rollout & 60%; 30 Jan 2025.
5. [Council — Jan 2025 public reports pack](https://edemocracy.coventry.gov.uk/documents/g13324/Public%20reports%20pack%2030th-Jan-2025%2010.00%20Communities%20and%20Neighbourhoods%20Scrutiny%20Board%204.pdf?T=10) — original council PDF; 30 Jan 2025.
6. [Council — FOI648678140 stock and hours](https://www.coventry.gov.uk/downloads/file/43696/foi648678140-redacted) — FOI original PDF; 8 Oct 2024.
7. [Council — new LED programme starts in 2025](https://edemocracy.coventry.gov.uk/mgAi.aspx?ID=44974) — primary 2025 decision; Oct 2025.
8. [Council — lighting programme overview](https://www.coventry.gov.uk/street-scene/streetlighting) — current status including Sept 2025 reversal; current.


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.
