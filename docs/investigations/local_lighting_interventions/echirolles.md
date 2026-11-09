# Stage A case dossier — Échirolles

**Original NUTS2016 container:** FRK24 — Isère

**Mechanism:** Early operational switch-off before extensive LED conversion

**Intervention period:** 1 Aug 2022 (selected areas); 1 Oct 2022 (most roads)

**Status:** Primary documentary evidence found, but the precise treated lamp-by-lamp historical polygon and suitable internal controls are **not yet certified**. No VIIRS or ISS local result has been measured.

## Scientific test

A decline in VIIRS night radiance after October 2022 should be present on switched-off road/park footprints and weaker or absent on genuinely unswitched comparators, **before** many LED replacements occurred.

## Independently grounded treatment footprint and comparison

Provisional treatment: Échirolles municipal road lighting other than explicitly excluded main arteries; separately code municipal parks/gardens, business areas, public-building lights, and the July/Aug 2022 intervention near part of rue/square Eugène-Sue and the Pathé cinema car park. Comparator candidates: **avenue du 8 Mai 1945**, **cours Jean-Jaurès**, **Honhoué–Grugliasco–Kimberley** main axis, reported by city as kept on. Do not treat whole roads as confirmed untreated lamp columns; test land-use and existing brightness imbalance, and avoid comparison across 2023 LED refurbishment.

## Dated documentary evidence and caveats

Contemporaneous city notices set municipal programmes in July–September 2022, with the October 2022 road event. Current council page shows post-2025 01:00–05:00 and must **not** be back-projected: during 2022–24 roadway interruption was generally 00:00–06:00 except main roads. LED share was about **5% in 2022** and **40% in 2024**. The 2023 municipal newsletter names **120 heads on central cours Jean-Jaurès replaced in summer 2023**, indicating this operational control route later received LED intervention. A Dec 2023 metropolitan presentation contains meter `EP CITE poste 079` time series indicating reduction beginning Oct 2022; its action table appears to label 2023 dates inconsistently with contemporaneous city notices and graph: **use city notices for intervention dates**.

## Records required to lock Stage A footprint

Request 2022 approved switching maps, circuit IDs, column coordinates and exception list; installation/retrofit log by lamp from 2022–24; monthly/daily electricity use for feeder `EP CITE poste 079` and comparable always-on circuits; documentation about garden/parking ownership and time schedule; exact road-segment lengths. Try treatment window Oct–Dec 2022 before 2023 modernisation; snow-screen 2021 baseline and use actual VIIRS overpass timing.

## Original primary reports and corroborating documents

1. [Municipality — July 2022 energy sobriety announcement](https://www.echirolles.fr/actualites/un-nouveau-plan-de-sobriete-energetique-pour-la-ville) — council announcement; Jul/Aug 2022.
2. [Municipality — Sept 2022 public meeting & road switch-off timetable](https://www.echirolles.fr/actualites/le-plan-de-sobriete-partage-avec-les-habitant-es) — council announcement; Sep 2022.
3. [Municipality — July 2022 rue Eugène-Sue / cinema car park](https://www.echirolles.fr/actualites/extinction-de-leclairage-public-dans-la-rue-eugene-sue-et-dans-le-grand-parking-du) — street-specific announcement; 29 Jul 2022.
4. [Municipality — lighting hours, exclusions and LED progress](https://www.echirolles.fr/vie-quotidienne/cadre-de-vie/extinction-de-leclairage-public) — current municipal status (2025 hours!); 2025/26 snapshot.
5. [Grenoble Métropole — Parole aux partenaires presentation](https://planclimat.grenoblealpesmetropole.fr/cms_viewFile.php?idtf=1784&path=Parole-aux-partenaires.pdf) — original slide PDF; includes meter graph; Dec 2023.
6. [Cité Échirolles no. 404 LED rollout](https://www.echirolles.fr/sites/default/files/2023-10/Cite%20Echirolles%20N404-web2_compresse.pdf) — city magazine PDF; 120 LED heads on Jean-Jaurès; Sep–Oct 2023.


## Analytical safeguards

- **NUTS2016 is a geographical container, not the treatment polygon.** Start from lamp/circuit/municipal documentation, never select treatments based on VIIRS brightening/dimming.
- Distinguish **street, column, lamp head, CMS node, feeder and VIIRS pixel**, including duplicate light heads and partial street coverage; prefer per-column timestamp and geometry.
- Call a source 'verified' only for what it actually establishes; current spatial inventories do not establish **historical** operation.
- Do **not** modify locked 2014–2024 European VIIRS, NUTS2016 cohort, 2017 zero baseline, 2021 snow correction, or fixed spider/moth PRC functions. Maintain standalone local outputs.
- Check actual nocturnal VIIRS overpass vs local switch-off clock time, BST/CET changes, snow/season, moon/cloud and effective observed-pixel coverage. Require pre-trends and compare within treated NUTS3 container, but examine spillover and treatment contamination. Existing 2020–22 monthly VIIRS collection is insufficient for several cases.
- ISS colour mosaics are irregular 2012–13 vs 2014–20 composites. No clean 2023–24 spectral inference from these composites; assess common-valid pixel counts and image dates before B/G and G/R comparisons.
- Energy use is not emitted luminous flux; B/G or G/R is not a unique spectral correction to DNB; no results are claimed before extracting local satellite pixels.
