# Stage A case dossier — Bristol

**Original NUTS2016 container:** UKK11 — Bristol, City of

**Mechanism:** Mixed LED + CMS upgrade, potential staggered/within-city analysis

**Intervention period:** Rolling citywide programme started March 2023; incomplete at Dec 2024

**Status:** Primary documentary evidence found, but the precise treated lamp-by-lamp historical polygon and suitable internal controls are **not yet certified**. No VIIRS or ISS local result has been measured.

## Scientific test

Compare streets with lantern+CMS replacement, previously LED streets receiving CMS-only nodes, and later/not-yet converted sites using actual installation dates; quantify the VIIRS signal of mixed spectral and operating change, not a pure LED coefficient.

## Independently grounded treatment footprint and comparison

Public Bristol streetlight spatial inventory has asset IDs, points/road and `CMS_CONTROLLED` present-day field. This is **not** historical switching record or retrofit dates. Define 2023/24 *historical* cohorts from City Council/Centregreat contract ledger: lantern replacement + node; node added to pre-existing LED; not-yet converted; heritage and excluded housing. Be careful: city council's project includes highways, some parks and harbour/car parks but excludes council housing. Within UKK11 match land-use and baseline brightness; do not use mixed broad city average as intervention.

## Dated documentary evidence and caveats

Council states rolling LED/CMS scheme began March 2023 across ~36,000 lights, aiming end 2025. Council's 2026–28 climate plan retrospectively reports **29,366 lanterns+nodes and 4,500 nodes-only** in completed programme; these totals describe work after 2024, **not** 2024 statuses. Supplier says 38 schedule variants and half-hourly energy/activity telemetry by asset (historical availability unverified). Original ISS composite ends by 2020 so not suitable for 2023–24 changes.

## Records required to lock Stage A footprint

Request council/contractor asset-ID geocoded historical snapshots (as at 2022, 2023 and 2024), old lamp and replacement type/CCT/SPD/watt/lumens, actual 2023–24 commissioning timestamps, `CMS_CONTROLLED` historical change, dimming schedules and 00:00–03:00 half-hour telemetry, maintenance outages. Council streetlight team: lighting@bristol.gov.uk (contact published in council design guidance). DNB 2022–24, respecting incomplete endpoint and per-pixel historical condition.

## Original primary reports and corroborating documents

1. [Bristol Council — LED upgrade project](https://www.bristol.gov.uk/residents/streets-travel/transport-plans-and-projects/upgrading-to-led-street-lighting) — primary timeline and scope; 2023–25.
2. [Bristol Council — Climate Emergency Action Plan 2026–28](https://www.bristol.gov.uk/files/documents/5241-climate-emergency-action-plan/file) — PDF; retrospective programme totals, p.11 and appendix; 2026.
3. [Schréder — Bristol bridging Smart City gap](https://uk.schreder.com/en/insights/blogs/bristol-bridging-smart-city-gap) — supplier primary account of CMS telemetry and schedules; 17 Apr 2026.
4. [Bristol Council — streetlight public GIS catalogue](https://opendata.bristol.gov.uk/datasets/bcc::streetlights-3) — asset points IDs + CMS flag; site intermittently inaccessible; current.
5. [Bristol Council — streetlight GIS REST layer](https://maps2.bristol.gov.uk/server2/rest/services/ext/ll_transport/MapServer/37) — GIS endpoint; direct access from web reader unavailable; current.
6. [Bristol Council — street lighting design guide/contact](https://services.bristol.gov.uk/files/documents/4350-03-8-2-street-lighting) — contact and specification PDF; current.


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.
