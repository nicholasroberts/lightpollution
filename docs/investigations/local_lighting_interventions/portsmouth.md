# Stage A case dossier — Portsmouth

**Original NUTS2016 container:** UKJ31 — Portsmouth

**Mechanism:** Sequential LED replacement then adaptive dimming

**Intervention period:** 2018–19 LED replacement; 2019 pilot; summer 2020 adaptive rollout

**Status:** Primary documentary evidence found, but the precise treated lamp-by-lamp historical polygon and suitable internal controls are **not yet certified**. No VIIRS or ISS local result has been measured.

## Scientific test

Look for a spectral shift in 2018–19 associated with LEDs, followed by an *additional* change in DNB radiance on streets whose output is lowered in 2020. Original 2019 pilot had an intentional unchanged comparison area.

## Independently grounded treatment footprint and comparison

2019 trial named **Baffins, Cosham, St Jude, Milton and Fratton**; news report says one trial neighbourhood was an unchanged placebo/control. **Its identity and precise pilot boundaries are not established**: request original trial design. After citywide adaptive rollout, all dimmable LED lamps become candidate treated; city centre and exceptions remain higher output. Parts of **Old Portsmouth** and **Southsea** had heritage lamps not replaced. These are contextual comparators rather than automatically exchangeable controls.

## Dated documentary evidence and caveats

Primary council guidance confirms LED upgrading since 2018, 2019 pilot, citywide summer 2020 adaptive lighting. Contemporary Feb 2020 report describes most LED lighting dimmed ~50% in 22:00–06:00 window (check whether 50% is electrical input, output command or luminous flux). ~15,000 LED heads reported installed Feb 2020. Existing post-2013 ISS mosaic may contain pre/post-2018 imagery but has **irregular dates and cannot be used as an exact 2020 event marker**.

## Records required to lock Stage A footprint

Request 2019 trial area polygons and name of placebo area; original pilot evaluations, CMS schedules/command %, dates by column, luminaire specification and deployment logs 2017–20, exceptions, and council geocoded lighting inventory. Prioritise 2018–21 monthly DNB, exact VIIRS overpass relative to 22:00–06:00 schedules.

## Original primary reports and corroborating documents

1. [Portsmouth Council — lighting replacement decision history](https://portsmouth.moderngov.co.uk/mgIssueHistoryHome.aspx?IId=11501) — council programme entry; 2017.
2. [Portsmouth Council — stray/unwanted streetlight guidance](https://www.portsmouth.gov.uk/wp-content/uploads/2022/06/Guidance-for-dealing-with-stray-and-unwanted-light-from-Street-lights-0622.pdf) — original PDF confirms two-stage timeline; Jun 2022.
3. [Portsmouth News — 2019 trial of dimming](https://www.portsmouth.co.uk/news/politics/street-lights-in-portsmouth-set-to-be-dimmed-in-bid-to-save-cash-1310884) — contemporaneous secondary report; 28 Feb 2019.
4. [Portsmouth News — five trial neighbourhoods, unchanged control and citywide rollout](https://www.portsmouth.co.uk/news/politics/this-is-why-street-lights-across-portsmouth-will-be-dimmed-at-night-1384884) — contemporaneous reporting; 5 Feb 2020.


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.
