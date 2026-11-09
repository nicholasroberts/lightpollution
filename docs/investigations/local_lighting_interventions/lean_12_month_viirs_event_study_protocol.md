# Lean 12-month VIIRS event-window protocol — five lighting interventions

**Prepared 9 Nov 2026; planning note, no local satellite extractions or effects yet.** Replaces the assumption that column-by-column histories are a prerequisite for a **first-pass municipality-scale VIIRS test**. Historical asset-level records remain **optional explanatory refinement**, not required input for exploratory descriptive trends. The NUTS2016 container and the locked 2014–24 European annual/PRC work remain unchanged.

## Resolution and inferential unit

- **VIIRS/SNPP DNB**: about **750 m native effective instrument resolution**. Earth Observation Group monthly V1 tiles have **15 arc-second** gridded cells (roughly **463 m north–south** and **~285–300 m east–west** at 50–52°N, corresponding to ~0.13–0.14 km² per grid cell). The grid spacing **is not independent optical resolution**: adjacent pixels correlate due to the native sensor footprint, geolocation, atmospheric scatter and processing.
- The five municipalities have **tens to hundreds of gridded cells**, not just one to a few, for their complete boundaries. But **independent scene-level information is smaller**, and narrow street-level changes may be diluted below detection. Exact count of total, lit and cf_cvg-valid pixels must be computed from raster and actual boundaries, not guessed from area.
- **ISS original astronaut photos** can be much finer spatially than VIIRS, but georeferenced archival 2012–13 / 2014–20 mosaic resolution and valid common-pixel coverage are **not yet inspected locally**. More importantly, the later mosaic is an **irregular collection of dated images**, not a monthly series. It cannot measure Échirolles 2022, Coventry 2024 or Bristol 2023–24 events. ISS is *optional spectral context*, not compulsory for the 12-month event test.
- **Time of acquisition**: Suomi NPP approximately 01:30 local *solar* time; actual cloud-free scenes and daylight saving may vary. Do not infer an actual switch-off for a site if the acquisition preceded the local switch-off. EOG monthly data pool eligible acquisitions and have no within-month intervention date information.

## Minimum independent documentary data needed now

For each case only:
1. Reliable **intervention start month / implementation window**; tag a month as transition/exclude if implementation extended across it.
2. Broad **geographic extent** (citywide, neighbourhoods, listed streets, municipal roads or affected lighting share) adequate to select a predeclared municipal/urban-core mask.
3. **Intervention mechanism and overnight operating hours**: switch-off; dimming; LED conversion; restoration of all-night lighting; combinations.
4. Evidence of **major competing interventions** during the same ±6-month window.

Do **not** wait for IDs/coordinates for every column, precise watts/lumens, CMS logs or per-day commissioning to make the first plots. Request those for mechanistic interpretation if the city-level effect is mixed/weak or if a meaningful treated-versus-control spatial contrast requires them.

## Primary temporal design

- Display a **12-month centred window**, nominally six pre- and six post-intervention calendar months. Mark *mixed/transition months separately*, and do not count these as unexposed or completely treated. In a staged rollout, plot rollout progress if independently available; an exact midpoint is a reference marker, not a real instantaneous treatment.
- A naive mean of **six preceding months** versus the **six following months** is **not** the main causal test: it changes season by design. For each displayed month compute an anomaly against **the same calendar month of the prior year(s)**, and a corresponding difference to a locally plausible comparator. Prefer at least one full pre-intervention annual cycle, or two where available. Use unchanged data sources and monthly definitions, with fixed geographic mask.
- Graph both *raw monthly radiance* and *seasonally referenced contrast*, and inspect pre-trends rather than automatically fitting an intervention model. A within-municipality comparand may be structurally different (major roads vs residential) and should be labelled as such; a neighbouring urban area with similar pre-trends may provide a second control when unaffected.
- **Coventry:** rollout occurred **May–July 2024**; mark all three transition months, use Jan–Apr 2024 pre and Aug–Dec 2024 after, with 2023 same-calendar-month reference. The 2014–24 frozen endpoint prevents a true symmetric six-month post-window if it would extend into 2025. Do **not** extend the primary study past Dec 2024 merely to achieve window symmetry. A separate 2025 sensitivity could be proposed later only as distinctly labelled, since Coventry lighting was restored in Sep 2025 and LEDs began Nov 2025.
- **Échirolles:** primary roads event **1 Oct 2022**, but LED conversion from 2023. Use April–Sept 2022 vs Oct 2022–Mar 2023 *as the 12-month display*, with transition flag where appropriate. Snow/ground reflectance, especially winter 2022–23, are crucial.
- **Portsmouth:** do *two* panels for 2018/19 LED rollout and **summer 2020** adaptive dimming. LED rollout is staged without one instantaneous date; treat it as gradual unless an exact phase start is documented.
- **Thanet:** first phase Dec 2016–May 2017, second phase from ~April 2018 subject to documented actual dates. Treat transitions as an **interval**, not a single abrupt breakpoint, and only hypothesise brightening until measured.
- **Bristol:** begins March 2023 and continues beyond 2024, so a ±6-month *start* window is primarily a **partial rollout** rather than a post-completion effect. Test 2023 and 2024 citywide trajectory, optionally predeclare geographic groups from the contemporaneous street or ward progress information. If not available, interpret only a descriptive trend during rollout.

## Monthly radiance quality rules

- Use the **original monthly EOG avg_rade9h and cf_cvg** inputs on a *fixed* raster grid and *fixed* comparison masks, not averages of annual NUTS3 parquet outputs. The pre-existing repository has 2020–22 both **VCMCFG/VCMSLCFG** monthly layer sets; additional targeted months needed for other events. Keep product configuration stable across before/after; compare the other product as a sensitivity analysis. **Do not mix annual average_masked v2 and monthly v1 as interchangeable values**.
- **cf_cvg is the number of acceptable observations contributing to a monthly pixel**, not a cloud fraction. A radiance zero with cf_cvg=0 is missing, not a physically dark observation. Show monthly median and interquartile range of cf_cvg, proportion of eligible pixels with cf_cvg>0, and sensitivity excluding low-coverage pixels. Test thresholds >=1, >=3, >=5 valid observations when enough support exists; do not set one arbitrary fixed threshold that silently drops winter months. Report days/actual covered scenes where product metadata allow.
- EOG masks clouds and certain moon/twilight/stray-light observations, but monthly coverage remains variable with weather and season. The **cloud-free observed nights may be systematically different nights** across months, e.g. weekends, temporary lights or unusual weather. Monthly averaging does **not** automatically eliminate seasonality or weather effects. The 2021 winter snow/albedo anomaly previously discovered in this project is an explicit warning; inspect snowfall and surface state and, where needed, compare the accepted native-pixel surface correction, not the frozen annual value.
- For each area show total static VIIRS pixels, baseline lit/eligible pixel count, number observed each month, median/mean radiance, within-city heterogeneity, cf_cvg distribution, matched-pixel trend vs fixed-support trend, and comparator outcome. Calculate uncertainty using month-level blocks or appropriate spatial autocorrelation rather than treating hundreds of adjacent grid cells as independent treatment replicates.
- Quantify **absolute radiance differences first**; percentages at low-light baseline pixels can be unstable. Maintain 2017 zero-baseline and 2021 annual snow treatment *locked*; do **not** subtract a -0.235 offset again from monthly V1 scenes or assume their zero-level processing is identical to annual V2.

## Interpretation / decision gate

**Stage B pilot first:** For each city, identify a reproducible municipal mask plus an independently specified urban-core/treated-area approximation. Extract the 12-month VIIRS sequence and 12 months of seasonally matched reference (where available). Flag low cf_cvg, transition periods, overpass timing and snow. Publish raw paired tables and plots; assess *whether an event-like change exists at the city scale*.

- **If an unambiguous, robust discontinuity appears** and coincides with a verified operational switch-off, detailed per-lamp logs may *not* be necessary for the modest conclusion **that the event contains genuine operation-driven changes**. Additional engineering records still strengthen validation.
- **If a weak/mixed/null trend appears**, avoid equating no detectable citywide change with no local intervention; request finer footprints or acknowledge dilution / insufficient observations.
- **If the claim requires isolating HPS→LED spectral effects from dimming**, precise schedules and lamp types or dated ISS/spectral imagery *do* become essential; do not claim a numerical partition based on DNB alone.
- Current broad streetlists are enough to start **Coventry**; policy exception streets can be exploratory controls for **Échirolles**. For Bristol and Portsmouth the temporal focus may be more informative than trying to label individual lamp pixels before the first plots.

## Key references / verified project inputs

- EOG VNL V1 monthly specifications and cloud-free observation layer: https://eogdata.mines.edu/products/vnl/
- NASA SNPP VIIRS 750-m DNB sensor: https://modaps.modaps.eosdis.nasa.gov/services/about/products/viirs-c2-nrt/VNP02DNB_nrt.html
- NASA typical local overpass: https://science.nasa.gov/earth/earth-observatory/all-clear-on-the-western-front-147094/
- ISS calibrated 2012–13 / heterogeneous 2014–20 mosaics: https://zenodo.org/records/7677478 ; scientific methods https://pmc.ncbi.nlm.nih.gov/articles/PMC9473566/
- Repo: `src/diagnose_fixed_nuts2016_monthly_product_atlas.py` (`avg_rade9h`, `cf_cvg`, 2020–22 dual products), `src/validate_iss_spectral_shift.py` (paired valid ISS ratios), `docs/results_iss_spectral_vs_locked_viirs_2014_2024.md` (interpretation safeguards).
- Previous comprehensive evidence and contact requests are indexed in [README](./README.md).

**Decision:** Start with the **minimum data needed at VIIRS scale**, then escalate to council asset-level history only when the signal, scope or interpretation requires it. No quantitative local VIIRS trend has yet been calculated in this stage.
