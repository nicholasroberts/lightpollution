# First-pass lighting-source reconstruction for the 18 original PRC sites (11 November 2026)

**Scientific status:** Initial **historical 2018/2019 evidence register**, not a measured spectral inventory, not a fitted LED fraction and **NOT a change to locked spider/moth PRC functions**. The original 18 user-verified records (field dates and coordinates) are at [`data/reference/prc_spider_18_sites_verified.csv`](../../data/reference/prc_spider_18_sites_verified.csv). Per-site editable classification table: [`prc_calibration_sites_2018_2019_firstpass_lighting_mix.csv`](./prc_calibration_sites_2018_2019_firstpass_lighting_mix.csv).

## The critical distinction

1. **Streetlighting technology within a municipality** can often be documented retrospectively; it is *not* the same as **the weighted spectrum of all light contributing to sky brightness at the precise calibration site**.
2. Public road lights can be high-pressure sodium (HPS), low-pressure sodium (LPS), mercury, ceramic metal-halide, LED or multiple technologies. A binary HPS/LED classification is **at best a limited idealisation**.
3. University and private/commercial illumination may be LED-heavy even if public roads remain sodium-heavy. Neither municipal inventory counts nor proximity to a bright city centre identifies the nighttime **weighted skyglow source mix**.
4. The **dark PRC limit** may converge under HPS and LED assumptions if artificial contributions vanish and natural polarized sky dominates; in a saturating biological response, the **very bright PRC asymptote can converge as well**. The strongest spectral sensitivity can be at intermediate radiances. Don't impose such convergence without deriving the model.
5. There are **NO sites in London, Paris or Amsterdam** in these particular 18 original field records; the brightest ones are around central Oxford and a Wolverhampton supermarket. Do not label them LED without dated evidence.
6. Existing spider PRC transfer is **17/18 sites** and moth is **16-site fit**; original Slabhouse was excluded from the fits for the documented anomalous EOG zero, but must remain in this provenance register. Original `opensquare.csv` and `oxfordbussinesspark.csv` historical monthly values are identical despite different positions. Do not let the duplicate irradiance proxy determine local lamp type.

## Historical evidence and what it actually establishes

### Oxfordshire: **public streetlighting sodium-leaning in May 2019, private mix unknown**

- The county's **17 March 2020 January business-monitoring report**, p.14–15, recorded **14,132 LED conversions among 59,631 highway streetlights (23.7%) as of 31 January 2020**, including **3,000 during the 2019/20 financial year**. Most county-maintained lamps were thus still **non-LED even by January 2020**. This is *not* a per-site count, and non-LED ≠ confirmed HPS. [Official 2020 council report PDF](https://mycouncil.oxfordshire.gov.uk/%28S%2803u25dvpydye1mi1jfpj1cnn%29%29/documents/s50188/PSC_MAR1220R05%20-%20Appendix%201%20-%20January%20Business%20Management%20Monitoring%20Report.pdf).
- Importantly, LEDs **were present before** the larger replacement phase: a dated [24 August 2017 county report](https://news.oxfordshire.gov.uk/north-oxford-carterton-bicester-abingdon-and-didcot-are-the-focus-for-continuing-streetlighting-improvement-work/) describes replacement of 3,200 lamps including North Oxford neighbourhoods. The main [accelerated programme began 2021](https://news.oxfordshire.gov.uk/money-saving-streetlighting-project-set-to-be-completed-next-year/).
- The eight Oxford sites were all measured **14 May 2019**: `s05` MINI Plant Oxford, `s15` Museum, `s16` Open Square, `s02` Oxford Business Park, `s04` St Anne's College, `s18` Worcester Lawn, `s01` Worcester Pitch, `s06` Zoology. County-road background is **sodium-leaning / non-LED likely**, but *local private/campus/industrial luminaires are unverified*. The factory and business park require particular care because ownership is not county highway lighting.
- **Confidence:** moderate for *majority non-LED county road stock*, low for *actual astronomical skyglow spectrum at individual sites*.

### Wolverhampton: **municipal sodium-leaning in April 2019; private light unknown**

- [City Council, 23 June 2020](https://www.wolverhampton.gov.uk/news/wolverhamptons-street-lights-get-smart): started large-scale replacement of **27,570 sodium streetlights** with LED **after** all four original April 2019 field observations; [17 March 2021 council report](https://www.wolverhampton.gov.uk/news/thousands-street-lights-transformed-smarter-greener-city) dates scheme start to **April 2020** and states **7,164** replaced by March 2021.
- Some smaller LED schemes predated this: [Council, 26 March 2018](https://www.wolverhampton.gov.uk/news/wolverhampton-shines-bright-lighting-city-project-reaches-first-milestone) documents LED decorative/feature lighting in the centre, not conversion of the entire road network.
- Sites: `s07` Church Hill, `s09` Ounsdale, `s08` Penn Cricket Club, `s03` supermarket Waitrose/Tesco. Treat *public roads* as sodium-leaning in 2019, *supermarket car-park/facade* lights as unverified private sources, and *cricket/sports lights* as on/off dependent and unknown.
- **Confidence:** moderate-high for municipal sodium prevalence, low for weighted skyglow at each site.

### Somerset/Wells: **already mixed in 2018; dark rural sites not technology controls**

- [Somerset County Council *2018* Transport Plan (draft, April 2018), §4.4.2](https://somerset.moderngov.co.uk/Data/SCC%20Scrutiny%20for%20Policies%20and%20Place%20Committee/201804241000/Agenda/Appendix%20A%202018%20LTP%20Implementation%20Plan%20Draft%20v5.pdf) reports **24,715 LED highway lamps** and **28,200 remaining to be converted** at the county level. Roughly 47% of these counted lamps had been converted, but there is **no evidence for the specific Wells roads, private light sources or spectra**.
- `s13` East Horrington (original multi-year monthly VIIRS mean ~0.50), `s17` Slabhouse (~0.50), and `s11` South Horrington (~2.89) were measured **21 May 2019**. At the two darker sites the natural polarized sky could dominate, making a local HPS/LED distinction of limited influence on PRC, *but that must be model-checked*.
- **Confidence:** moderate for county-wide mix, low for local skyglow; no justified per-site HPS/LED percentages.

### Southern Sweden: **historical mix unknown**

- Simrishamn Municipality's [2016 public lighting maintenance tender](https://se.openprocurements.com/tender/2016-tjansteentreprenad-avseende-skotsel-och-underhall-av-offentlig-belysning/) names **both Borrby and Skillinge** as serviced places, and permits lamp replacements. It provides **no 2018 inventory of HPS versus LED**. Current municipal lighting pages are not an accurate substitute.
- `s14` Borrby Strandbad, 21 Sep 2018, low VIIRS radiance (~0.36); `s10` Skillinge, 21 Aug 2018, moderate (~5.56). Both **unknown technology**, though the Borrby site's low light level may give convergence in PRC sensitivity scenarios.
- **Confidence:** low. Do not assign sodium or LED dominance by national stereotype.

### Gavarnie: **dark site; historical lamp type unknown**

- `s12` Gavarnie (27 Sep 2018) was measured in a very dark Pyrenean setting, original multi-year mean ~0.56. The [Pic du Midi regional dark-sky reserve was recognised in 2013](https://www.cerema.fr/index.php/fr/presse/dossier/suivi-national-extinctions-nocturnes-cartographie-visualiser), but this establishes **dark-sky policy context, not local lamp emission spectra**.
- [Gavarnie-Gèdre municipal council discussed nighttime switch-off in October 2022](https://politique.pappers.fr/commune/document/e85da724a48ba709f5850dc924e8d283); **do not project that later policy back into September 2018**.
- **Confidence:** high for dark-context classification, low for technology.

## Per-site calibrated relationship (source radiance values are historic *original multi-year monthly means*, NOT locked observation-year calibrated annual VIIRS)

| Field site (ID) | Measurement date | Original VIIRS monthly mean¹ | Spider PRC | 2018/2019 *first-pass* source assignment |
|---|---|---:|---:|---|
| s18 Worcester College Lawn, Oxford | 2019-05-14 | 58.73 | 0.123 | Public sodium-leaning; college unknown |
| s01 Worcester College Pitch, Oxford | 2019-05-14 | 48.55 | 0.144 | Public sodium-leaning; college/pitch unknown |
| s16 Open Square, Oxford | 2019-05-14 | 41.83 | 0.149 | Public sodium-leaning; square sources unknown |
| s02 Oxford Business Park | 2019-05-14 | 41.83 | 0.192 | Public sodium-leaning; private park unknown |
| s15 Museum, Oxford | 2019-05-14 | 37.87 | 0.193 | Public sodium-leaning; museum/building unknown |
| s03 Waitrose/Tesco, Wolverhampton | 2019-04-21 | 35.36 | 0.122 | Public sodium-leaning; supermarket unknown |
| s04 St Anne's College, Oxford | 2019-05-14 | 29.04 | 0.157 | Public sodium-leaning; college unknown |
| s05 MINI Plant Oxford | 2019-05-14 | 25.20 | 0.165 | Public sodium-leaning; factory unknown |
| s06 Zoology, Oxford | 2019-05-14 | 19.61 | 0.236 | Public sodium-leaning; university unknown |
| s07 Church Hill, Wolverhampton | 2019-04-20 | 13.60 | 0.274 | Municipal sodium-leaning; other unknown |
| s08 Penn Cricket Club, Wolverhampton | 2019-04-22 | 10.10 | 0.262 | Municipal sodium-leaning; sports unknown |
| s09 Ounsdale, Wolverhampton | 2019-04-20 | 7.12 | 0.321 | Municipal sodium-leaning; residential unknown |
| s10 Skillinge, Sweden | 2018-08-21 | 5.56 | 0.276 | Technology unknown |
| s11 South Horrington, Wells | 2019-05-21 | 2.89 | 0.423 | County mixed; local unknown |
| s12 Gavarnie, France | 2018-09-27 | 0.56 | 0.470 | Dark setting; lamp type unknown |
| s13 East Horrington, Wells | 2019-05-21 | 0.50 | 0.419 | Dark/low artificial; county mixed |
| s17 Slabhouse, Wells | 2019-05-21 | 0.50 | 0.455 | Dark/low artificial; county mixed; excluded from fit |
| s14 Borrby Strandbad, Sweden | 2018-09-21 | 0.36 | 0.439 | Dark/low artificial; lamp type unknown |

¹ Native source `original_multiyear_monthly_radiance_mean` in `nW cm⁻² sr⁻¹`. It is **not** the independent 2018/19 EOG V2 annual field-year radiance used in the locked 17-site spider / 16-site moth fits. All locations, dates and PRC values from the source calibration reference; no new photometry, satellite values or manual LED fractions inferred.

## Implication for the next spectral response experiment

- The **bright end of the calibration** should not be called '100% LED', especially not given the explicit 2019 Oxford and Wolverhampton municipal evidence. It is better described as **an unknown, site-varying lighting mixture embedded in the existing empirically fitted curve**.
- The **dark end** is useful for constraining the natural-sky asymptote in a mechanistic model: if artificial skyglow is negligible, the predicted cue contrast becomes insensitive to the *artificial lamp spectrum*. The original empirical spider and moth models may already nearly saturate at **very bright radiances**, which could make spectral uncertainty greatest in the middle.
- The site-by-site private/muncipal mixture must not be inferred from a county-wide percentage. To identify lamp technologies locally, next evidence should be (a) time-specific university/supermarket industrial luminaire inventories or streetlight GIS lamp-type attributes, (b) 2018/19 ISS RGB imagery near sites, and (c) original experiment logbooks/night-time field photos or spectra if available.
- **Recommendation:** first model spectral uncertainty at all sites using *plausible* lighting scenarios; avoid locking the empirical curve to an assumed HPS or LED class without direct independent evidence.
