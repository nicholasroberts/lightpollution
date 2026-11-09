# Stage A: independent treatment footprints and comparators

**Prepared 9 November 2026.** No local effect estimation yet. All IDs/treatments below must be assigned from lighting records **before** extracting VIIRS outcomes. Confidence legend: **documented policy**, **published street list**, **published lamp points (current only)**, **unresolved per-lamp histories**.

| Locality | Treatment footprint that can now be nominated independently | Candidate comparison | Independently available spatial/operating source | Essential missing fact before locked analysis | Stage A status |
|---|---|---|---|---|---|
| Échirolles | Most municipal roads under Oct 2022 midnight–06:00 switch-off; separately Aug parks/business sites and rue Eugène-Sue/Pathé car park | City-named always-on axes: avenue du 8 Mai 1945; cours Jean-Jaurès; Honhoué–Grugliasco–Kimberley | 2022 notices and 2023 independent municipal feeder energy time series (`EP CITE poste 079`) | 2022 streetlight circuit/lamp GIS, exact switching and LED conversion dates; in 2023 Jean-Jaurès gained 120 LEDs | **Named road exception footprint**, no authoritative point geometry |
| Coventry | Official council 51-page switched-street + ward list, mirrored as source PDF/XLSX links; 2024 May–Jul switch-off | Verified 2024 exempt streetlight IDs, perhaps city-centre/priority roads; not simply all unlisted streets | Official street PDF/XLSX and scrutiny/FOI reports | Column geocodes, whether a whole street or only part switched off, exact May–Jul 2024 day, lamp technology and exemptions | **Published affected street list; no verified lamp-level geometry** |
| Portsmouth | 2019 pilot named Baffins, Cosham, St Jude, Milton, Fratton, later dimmed CMS areas | Unchanged pilot/placebo site (unidentified); exempt/heritage lamps with caveats | Council 2022 PDF timeline and Feb 2020 contemporaneous trial coverage account | Identify which pilot area was unchanged; street/column pilot map; dimming % as luminous or electrical output; column-level CMS activation | **Five named pilot neighbourhoods; control identity missing** |
| Thanet | Kent Dec 2016–May 2017 converted stock vs 2018 delayed concrete-column stock, plus all-night restored operation | Later-converted stock or separately documented always-on/part-night lamps, subject to baseline matching | County conversion policy, historical counts, **current Kent lamp location CSV and ArcGIS points** | 2015–19 ID-dated lamp technology + actual all-night commissioning dates and prior switch-off schedule | **Current points available; history missing** |
| Bristol | 2023–24 lantern+CMS replacements; existing LED nodes-only; not-yet-converted sites | Not-yet-converted before their actual date; nodes-only as mechanistically distinct treatment | Public Bristol streetlight GIS with IDs/current CMS flag, council 2023 plan, 2026 retrospect & CMS supplier report | Asset ID–linked 2023–24 actual swap/commissioning dates, previous technology, light output, dimming schedules and energy telemetry | **Current asset points available; historical treatment dates missing** |

## Classification schema (future CSV/GeoPackage)

`case_id, nuts2016_id, boundary_version, source_id, asset_id, lamp_id, circuit_id, geometry_source, lamp_position, owner, old_technology, new_technology, prior_cct, post_cct, prior_watts, post_watts, prior_lumens, post_lumens, event_type, scheduled_start_local, actual_commissioned_local, before_schedule, after_schedule, overpass_status, treatment_group, control_group, confidence_grade, source_url, source_accessed_date, notes`.

## Spatial rule

Build *lamp point* and street segment joins where records allow; do not equate the full municipality to a treatment polygon. Generate exposure intensity for a native VIIRS pixel as a function of recorded columns/lamps and influence radius, report sensitivity to buffer 0/250/500 m and to control contamination, and preregister before seeing local DNB trends. These buffer sizes are proposed **sensitivity settings**, not calibrated physical light-spread radii.

## Prerequisite acceptance tests

1. An independent documentary record supports treatment type and **date** with a locality, street, circuit or asset identifier.
2. Actual physical geometry or location-to-centreline geocoding has been reviewed and NUTS2016 geography verified; older/current administrative boundaries not conflated.
3. Comparison is plausibly unaffected at the test date (or explicitly documented as different intensity) and pre-trend inspectable.
4. Actual VIIRS overpass falls inside the operational change window; annual composites' calendar mixing clearly reported.
5. City-centre/major-axis/heritage controls are flagged as *structurally different* and used as exploratory controls unless more comparable untouched streets are found.
6. Public current GIS attribute `CMS_CONTROLLED` or any current LED status **never retroactively labels a 2016/2024 lamp**.

## Source links

See [the full primary-document catalogue](./primary_source_reports.md) and five [case dossiers](./README.md).


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.
