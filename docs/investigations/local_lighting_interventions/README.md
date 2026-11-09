# Local lighting interventions — five-case evidence hub

**Prepared:** 9 November 2026 · **Branch:** `publication-figures` · **Status:** Stage A documentary footprint investigation; NO local VIIRS/ISS effect estimates yet.

This folder consolidates our original five-case register, independently authored research briefs and source-report catalogue. The study tests whether changes in DNB radiance can arise from independently documented operational switching/dimming, spectral LED conversion and/or longer lighting hours. It is supplementary validation of the locked European study; it does **not** refit PRC curves or select treatments by satellite outcome.

## Comparative motivation

| Case | NUTS2016 | Intervention during locked 2014–24 window | Mechanistic comparison | Main strength | Main limitation |
|---|---|---|---|---|---|
| **Échirolles** | FRK24 Isère | Overnight switch-off Oct 2022; significant LED rollout after 2022 | Primarily operating hours **before** LED change | City-named permanently illuminated exceptions; independent municipal feeder electricity graph | Streetlight circuits and exact lamp points not yet geocoded; 2023 changes contaminate later controls |
| **Coventry** | UKG33 | Part-night switch-off implemented May–Jul 2024; main LEDs only Nov 2025 onward | Switch-off without contemporaneous *principal* LED conversion | **51-page switched-street list + XLSX** and council timeline | Lamp-by-lamp exemptions/dates and 2024 LED mix missing; city centre controls not exchangeable |
| **Portsmouth** | UKJ31 | 2018–19 LEDs, 2019 dimming pilot, citywide adaptive dimming summer 2020 | Spectral conversion then subsequent output change | Named trial neighbourhoods with reported unchanged placebo area | Trial placebo identity and dimming logs not yet obtained |
| **Thanet** | UKJ44 East Kent | LEDs from late 2016, move to optimised all-night operation | Potential **brightening** from longer hours despite LED conversion | Early/delayed cohorts; county current lamp inventory + GIS | Exact restoration/commissioning dates absent; delayed columns not random |
| **Bristol** | UKK11 | Rolling LED + CMS started Mar 2023, ongoing after Dec 2024 | Lantern+node vs node-only changes, staggered cohorts | Public geocoded assets, possible half-hour CMS telemetry | Historical 2023–24 schedules and install logs not publicly available |

## Documents here

- [Original five-case research register](./evidence_register_2014_2024.md) — consolidated context, original sources, satellite windows, limitations.
- [Stage A footprint matrix and decision rules](./stage_A_footprints.md) — what treatment and controls can actually be defined, and what is missing.
- [Primary source reports and original-document catalogue](./primary_source_reports.md) — one index of source reports/PDFs/XLSX and council records, with original URLs. **We link original reports, rather than reposting third-party PDF/Excel documents without verified redistribution rights.**
- [Data requests to local authorities](./data_requests.md) — targeted requests for asset/circuit IDs, commissioning and CMS logs.
- Five evidence dossiers: [Échirolles](./echirolles.md) · [Coventry](./coventry.md) · [Portsmouth](./portsmouth.md) · [Thanet](./thanet.md) · [Bristol](./bristol.md).

## Status of source documents

This folder contains **our complete analytical research documents** in GitHub. The external council, supplier and news report *originals* remain at source URLs and are catalogued with their formats, dates and uses. They have **not** been copied as PDFs into this repository: redistribution licences were not established, and the GitHub text-file workflow does not preserve those binaries. Do not mistake links for locally archived reports.

## Further reading

- Previously established ISS datasets: https://zenodo.org/records/7677478
- Existing spectral pipeline: `src/validate_iss_spectral_shift.py`, `src/analyze_iss_spectral_vs_locked_viirs.py`, `docs/results_iss_spectral_vs_locked_viirs_2014_2024.md`


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.
