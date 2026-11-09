# Stage A case dossier — Thanet

**Original NUTS2016 container:** UKJ44 — East Kent

**Mechanism:** Opposite-direction hypothesis: LED retrofit plus restoration of all-night operation

**Intervention period:** County policy Mar 2016; Thanet conversions Dec 2016–May 2017, further phase from 2018

**Status:** Primary documentary evidence found, but the precise treated lamp-by-lamp historical polygon and suitable internal controls are **not yet certified**. No VIIRS or ISS local result has been measured.

## Scientific test

Find whether reinstating night operation on previously part-night streets produces increased DNB radiance, notwithstanding VIIRS spectral sensitivity to LEDs. Brightening is **hypothesis only**, not measured result.

## Independently grounded treatment footprint and comparison

Thanet district (Margate, Ramsgate, Broadstairs etc.) within UKJ44. Kent's public streetlight location CSV includes Thanet; a county GIS MapServer layer contains points and asset IDs. Neither current dataset confirms historic lighting technology or actual switching time. 2017 inventory reported **9,931** lamps, **6,922** converted in Dec 2016–May 2017 and **3,009** outstanding; around **2,500** required column changes. Request historic per-column phase/geocodes and old part-night vs all-night schedule. Within Thanet, compare early 2016–17 retrofit with delayed concrete-column stock, but 'delayed' may differ in road type or prior operating hours.

## Dated documentary evidence and caveats

County policy approved March 2016 allowed optimised all-night as new LEDs entered CMS. A 2014 committee record shows preceding part-night policies had involved switch-offs between around midnight and early morning. Do not assume identical operation for every Thanet street; seek original commissioned schedule. Kent Sept 2018 committee minutes document concrete-column delays; a planned April 2018 phase 2 is **not proof all installations completed then**. Some assets may belong to other authorities. Geocoded present-day points must not be assumed historical.

## Records required to lock Stage A footprint

Retrieve Kent published Thanet-inclusive location CSV; request 2015–19 inventory snapshots and ID–date–lamp type–CMS start logs, records of exact all-night restoration and previous part-night hours, and concrete-column replacement maps. Compare 2014–15 vs 2016–19 monthly DNB by declared cohorts with coverage corrections, controlling spill from adjacent Kent schemes.

## Original primary reports and corroborating documents

1. [Kent County Council — 2016 policy decision](https://democracy.kent.gov.uk/ieDecisionDetails.aspx?Id=878) — policy record with source reports and appendices; 4/15 Mar 2016.
2. [Kent Council — 2014 part-night programme minutes](https://democracy.kent.gov.uk/ieListDocuments.aspx?CId=831&MID=5642&txtonly=1) — historical policy minutes; 24 Apr 2014.
3. [Thanet News — counts converted / outstanding](https://theisleofthanetnews.com/2017/07/19/kcc-has-been-unable-to-convert-a-staggering-2500-streetlights-in-thanet-to-led-all-night-lighting/) — reported council figures; 19 Jul 2017.
4. [Thanet News — planned second phase](https://theisleofthanetnews.com/2018/03/05/conversion-of-3009-thanet-streetlamps-to-led-all-night-lighting-to-begin-in-april/) — 2018 reporting; projected schedule only; 5 Mar 2018.
5. [Kent County Council — Sept 2018 committee and minutes](https://democracy.kent.gov.uk/ieListDocuments.aspx?CId=831&MID=7976) — primary concrete-column delay confirmation; 20 Sep 2018.
6. [Kent County Council — streetlight location CSV for Thanet](https://www.kent.gov.uk/roads-and-travel/what-we-look-after/street-lights/faulty-street-lights) — public asset IDs and location CSV; current snapshot.
7. [Kent Council ArcGIS — street lights point layer](https://webapps.kent.gov.uk/server/rest/services/Test/KentViewPortalTest/MapServer/44) — asset GIS spatial layer; temporal metadata absent; current snapshot.
8. [Kent County Council — Highways Asset Management Plan](https://www.kent.gov.uk/__data/assets/pdf_file/0003/124266/Highways-Asset-Management-Plan.pdf) — post-intervention contextual PDF; 2024.


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.
