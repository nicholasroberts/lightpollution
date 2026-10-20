# Evidence review: why VIIRS night-light radiance declined in parts of Europe (2014–2024)

**Research note — 20 October 2026. Status: literature/context review, NOT causal attribution of the project's NUTS3 radiance change.**

## Research question and boundary

The project's corrected 2014–2024 VIIRS radiance maps show particularly extensive declines in parts of France, the UK, Belgium and the Netherlands, with contrasts against Germany, Switzerland and Italy. A decline in the **VIIRS DNB band** is not, on its own, evidence that broad-spectrum artificial-light emission, biological exposure or skyglow fell proportionally.

Three non-mutually exclusive mechanisms need distinguishing:

1. **Real operational changes:** lamps switched off at the satellite overpass, part-night lighting, dimming, lower luminous flux, changed spatial shielding/directionality, reduction of installed lighting.
2. **Spectral substitution:** high-pressure sodium → white LED changes how much emitted light falls into the DNB's relatively insensitive short-wavelength/blue part of the spectrum. A decline in DNB radiance can occur even when the visible/scotopic/biologically active light does not fall.
3. **Environmental/measurement effects:** cloud, atmosphere, humidity, viewing geometry, snow/surface albedo, stray-light/background offset, annual processing and mask changes, and changing satellite/source products. Our **2017+ 0.235 baseline offset** and spatial Jan–Feb **2021 snow/albedo candidate** tackle specific parts of this third mechanism, NOT the others.

**Warning about independence:** Cerema's French switch-off inventory uses VIIRS, so it is **independent institutional analysis but not independent sensor corroboration**. The April 2026 *Nature* paper uses daily VIIRS-derived NASA Black Marble, independent processing/temporal method but a **shared spectral sensor**. The Enedis overnight electricity/power observation and legal/municipal policy records are genuinely different evidence types.

## France — unusually strong real-world policy evidence

**Primary government source:** Cerema, *Extinction de l'éclairage public : une cartographie sur les pratiques des collectivités* (30 July 2025), https://www.cerema.fr/fr/actualites/extinction-eclairage-public-cartographie-pratiques and formal dossier https://www.cerema.fr/fr/presse/dossier/suivi-national-extinctions-nocturnes-cartographie-visualiser .

Cerema's retrospective **2014–2024** VIIRS analysis covered **19,262 analysable French municipalities**. Within that sample, **11,980 (62%)** were classified as practising complete middle-of-night switch-off; **3,547 (18%)** partial switch-off **or major refurbishment**. These are percentages **of the analysed sample**, NOT of all municipalities or of emitted radiance, and the combined partial/refurbishment classification is not a proof of extinguishing the corresponding lamps. Cerema estimates **30% of total switch-offs were implemented during September–December 2022** energy-crisis period. Small municipalities were more likely to have complete switch-offs; larger ones tended to have partial/refurbishment patterns. A subsequent 2026 update for **2024–2025** reported 73% of 22,773 examined municipalities engaged in some kind of full/partial switch-off: https://www.cerema.fr/fr/actualites/extinction-eclairage-public-nouvelles-donnees-2024-2025 . The 2025 figure is **outside our 2024 endpoint** and must not be used as the 2014–24 statistic.

**Independent electricity-side observation:** Cerema summarises Enedis data showing a **34.2% drop in average street-lighting electric power drawn around 02:00**, between **May 2022 and May 2023** (not a 34.2% cut in full-year electricity use or optical radiance): https://www.cerema.fr/fr/activites/mobilites/chiffres-cles (September 2023 item). Supports operational night-time changes, but the measurement's sectoral/spatial support differs from our NUTS3 satellite means.

**Legal source:** France's **27 December 2018 arrêté** imposed restrictions for multiple types of commercial, non-residential, shopwindow, heritage and park lighting and constrained upward emission and colour temperature: https://www.legifrance.gouv.fr/loda/id/JORFTEXT000037864346/ , explained by Ministry of Ecology https://www.ecologie.gouv.fr/politiques-publiques/arrete-du-27-decembre-2018-relatif-prevention-reduction-limitation-nuisances . Do **not** call this a nationwide mandatory midnight switch-off of all streetlights: widespread municipal streetlight switching is a distinct phenomenon.

**Interpretation:** Strong support for *some genuine operational dimming/switch-off* contributing to the VIIRS decline, in addition to LED spectra and observational effects. Attribution of our exact NUTS3 effect sizes requires matching implementation dates/municipal boundaries and independent electricity/lighting inventories, not citing national percentages as region-level correction multipliers.

## United Kingdom — part-night lighting and LEDs act together

**Peer-reviewed empirical history:** Steinbach et al. (2015), *The effect of reduced street lighting on road casualties and crime in England and Wales*, Journal of Epidemiology and Community Health; https://doi.org/10.1136/jech-2015-206012 . Of **62 English/Welsh participating local authorities with usable data** (not a national random percentage), **48%** had introduced part-night lighting, **65%** dimming, **84%** white-light conversion and **8%** some switch-off by the analysed 2000–2013 period. These strategies started substantially before 2014.

**Direct local examples (establishment of mechanisms, NOT national attribution):**
- Hampshire County Council: documented dimming since ~2010, part-night lighting introduced 2019: https://www.hants.gov.uk/aboutthecouncil/haveyoursay/consultations/future-services-consultation/street-lighting .
- Coventry Council: part-night lighting programmed in 2024 for ~60% of eligible streetlights: https://edemocracy.coventry.gov.uk/mgAi.aspx?ID=43074 .
- Bristol City Council: ~36,000 streetlights in LED conversion programme begun March 2023, with programmable dimming; https://www.bristol.gov.uk/residents/streets-travel/transport-plans-and-projects/upgrading-to-led-street-lighting .
- Lincolnshire Council: LED, dimming and part-night conversion project begun 2016; https://www.lincolnshire.gov.uk/news/article/1669/street-lights-by-the-numbers .

**Independent chromatic evidence:** Sánchez de Miguel et al. (2022) reports strong shifts in UK blue/green and/or green/red ISS ratios; approximately 51% of public UK streetlights had reportedly converted to LEDs by early 2019 according to cited national sources. See spectral section below. LED conversion and actual switching/dimming could jointly drive the UK's measured trend; do not assign fixed fractions without more evidence.

## Germany — a comparison, not a zero-change control

The September 2022 German energy-saving ordinance and amendments constrained illuminated advertising **22:00–06:00** and some public-building/monument lighting: German Federal Economics Ministry, https://www.bundeswirtschaftsministerium.de/Redaktion/DE/Pressemitteilungen/2022/09/20220928-kabinett-beschliesst-anpassungen-bei-energieeinspar-verordnung.html ; original temporary regulation with scheduled expiry 28 February 2023 at https://www.bgbl.de/xaver/bgbl/start.xav?start=%2F%2F%2A%5B%40attr_id%3D%27bgbl122s1446.pdf%27%5D . These temporary rules cannot explain a smooth decade-long 2014–24 trend by themselves.

Sánchez de Miguel et al. (2022) found **less pronounced ISS-colour spectral changes in Germany and Austria** than in the UK/Italy/Ireland/Romania. That makes German regions useful **contrasts in LED timing and policy**, but Germany must not be labelled an unchanged or untreated control: municipal LED modernisation and local restrictions varied.

## Switzerland — official satellite indicator and local-policy variation

Swiss Federal Office for the Environment national indicator *Light emissions* (updated April 2026): https://www.indikatoren.admin.ch/public/v2/detail?ind=NI004&lng=en . VIIRS-based night-light indicator **increased by 7.91%** comparing mean 2014–16 with mean 2018–20. This metric differs in years, footprint and processing from our fixed NUTS2016 2014/24 pair, and shares VIIRS spectral sensitivity; it does **not** independently prove the underlying emitted visible light increased by 7.91%. Swiss local practices can differ materially; avoid a blanket "Swiss non-dimming control" claim.

## Italy — strong spectral-change caution

Kyba et al. (2017), *Artificially lit surface of Earth at night increasing in radiance and extent*, *Science Advances* 3:e1701528, https://doi.org/10.1126/sciadv.1701528 , full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC5699900/ . Figure 5 documents **Milan's sodium-to-white-LED change between ISS images in 2012 and 2015**, accompanied by **lower VIIRS DNB radiance despite the visible spectral change**, because DNB largely lacks 400–500 nm sensitivity. The authors warn that ~30% VIIRS declines *could* be caused by complete HPS→LED conversion without proportional loss of visible light; this is an illustrative technology scenario, **not a conversion factor to apply to all regions**.

Italy is therefore especially appropriate for identifying how spectral change and VIIRS radiance can diverge; the timing of Milan's change partially predates our 2014 baseline, so do not assign all Italian NUTS3 2014–24 trends to this specific conversion.

## Cross-country independent processing: 2026 Nature study

Li, T. et al. (2026), *Satellite imagery reveals increasing volatility in human night-time activity*, *Nature* 652, 379–386, https://www.nature.com/articles/s41586-026-10260-w , published 8 April 2026. Analyses daily **NASA VIIRS Black Marble** with change detection, **2014–2022** (NOT 2014–24), finds net radiance reductions **France ~33%, UK ~22%, Netherlands ~21%** within their study definitions, and **~4% net dimming in Europe**. It explicitly discusses LED conversions, lighting restrictions and energy conservation. The time-series method supports widespread true changes in the *sensor record* and regional heterogeneity. It **does not overcome the common DNB spectral blind spot**, nor isolate physical light-output causes at our NUTS3 scale. Do not numerically equate these relative change figures to ours without harmonising source product, denominator, geography and years.

## ISS spectral ratios: main next methodological advance

Sánchez de Miguel, A., Bennie, J., Rosenfeld, E., Dzurjak, S. & Gaston, K. J. (2022), *Environmental risks from artificial nighttime lighting widespread and increasing across Europe*, *Science Advances* 8:eabl6891. DOI https://doi.org/10.1126/sciadv.abl6891 ; accessible full text PDF https://zenodo.org/record/8194184/files/sciadv.abl6891%20%282%29.pdf .

The study uses colour-calibrated astronaut ISS photographs to compare broad **2012–2013 vs 2014–2020** epochs. Increases in **B/G and/or G/R** correspond to transitions away from sodium toward white LED lighting (and sometimes other white technologies); largest cited shifts include **Italy, Romania, Ireland and the UK**, with smaller shifts in **Austria and Germany**. Primary data openly available: https://zenodo.org/records/7677478 (version 1.02), https://zenodo.org/records/6468940 (large original epoch mosaics), processing code https://zenodo.org/record/4486581 .

**Critical limitations for our planned 'spectral correction':**
- The ISS colour comparisons are **two broad, irregular epochs**, not an annual 2014–24 spectral time series. They cannot supply exact NUTS3-by-year 2014–2024 LED correction factors without temporal assumptions/external conversion dates.
- Two broadband colour ratios alone generally do **not uniquely reconstruct spectral power distribution (SPD)**, luminaire luminous flux, radiated upward spectrum, skyglow or animal receptor stimulation. Need plausible HPS/LED SPD basis sets and published/actual instrument response functions.
- DNB roughly samples 500–900 nm; an ISS-inferred rise in blue light can lower DNB-to-full-visible conversion efficiency even as **scotopic** flux / short-wavelength scattering increases. A DNB-equivalent trend, a photopic/scotopic optical-flux trend, and an ecological PRC trend are **different estimands**; don't label one "true emitted light" without its spectral weighting.
- ISS colour data may be very sensitive to clouds, camera calibration, spatial resolution, illumination technology, seasonal/observation selection. Aggregate the high-resolution ISS layers with sensible masks to **fixed NUTS2016 NUTS3**, and propagate uncertainty.
- In the existing empirical PRC calibration the DNB R-to-contrast relationship is inferred from field visual-response data in mixed local lighting spectra. A new LED-corrected 'broadband R' predictor **cannot automatically be inserted into that pre-fitted DNB→PRC function** without checking whether the correction is defined on the scale that function represents. Compare baseline DNB-derived PRC with explicit spectral-response scenario(s), label the extra assumptions, and avoid double correction.

### Proposed bounded study: evidence chain not circular curve-fitting

A. Assemble per-NUTS3 2014–2024 **baseline-offset** and candidate **snow-adjusted** DNB annual trends; report the 2021 correction separately until signed off. No curve re-optimisation.

B. Map ISS colour ratios and **change in ratios** within the common European NUTS3 2016 geography; record data years, missingness and country. Compare with contemporaneous DNB trend directions and documented public-lighting conversion programmes.

C. Forward-model plausible lighting-source SPD mixes through (i) published VIIRS DNB response; (ii) ISS RGB channels; (iii) optionally photopic/scotopic and relevant biological spectral receptor efficiency. Infer a **bounded range** of physically possible correction to VIIRS detection of LED-associated blue emission rather than treating B/G or G/R as an exact correction coefficient. Avoid using a spectral multiplier to force monotonic brightness.

D. Re-evaluate **two separate empirical polarization-receptor-contrast species relationships (spider and moth)** on consistently defined observed DNB radiance and, as a separate sensitivity, on a physically defined spectral-light quantity if transfer comparability is justified. Pixel transform before NUTS aggregation. Do not import locked spider coefficients into moths.

E. Ground selected case-study NUTS3 regions (e.g., France vs German/Austrian examples, an Italian early LED transition, UK policy change) with **independent lighting schedules and electricity records**. In any country-year comparison, distinguish VIIRS-only policy inventory (not independent spectral evidence) from primary policy/electricity information.

F. Test population-density stratification to identify the **density/brightness ranges in which predicted spider and moth PRC changes are largest**, using baseline density measured independently of VIIRS wherever possible. Keep absolute biological ΔPRC, percentage ΔPRC and R trend distinct.

## Status / evidence hierarchy

| Mechanism | Evidence for existence | Can quantify fraction of OUR 2014–24 NUTS3 ΔR? |
|---|---|---|
| French overnight switch-off | Strong: Cerema classified municipalities; law; independent Enedis power reduction | Not yet: need geography/time matched |
| UK part-night/dimming/LED | Strong locally, policy-specific and peer-reviewed; ISS colour | Not yet: heterogeneous authority/year coverage |
| German national 2022 restrictions | Strong for narrow temporary time window | Not 2014–24 trend overall |
| Swiss radiance increase | Official VIIRS-derived national indicator | Different era and non-independent sensor |
| Italian LED spectral shift | Strong ISS/Milan + published DNB sensor demonstration | Not a uniform country-wide factor |
| Cross-European historical dimming | Strong satellite observation, corroborated by *Nature* Black Marble | Not yet equal to emitted broadband light |
| LED effect on species PRC | Plausible mechanism; test via paired spectroscopy/ISS + sensitivity | No direct validated ecological correction yet |

**Manuscript positioning:** the novel claim is not discovery that European VIIRS night lights can dim (documented elsewhere). It is (1) a transparent correction/interpretation of a decade of European night radiance including DNB offset, 2021 snow and LED spectral sensitivity, and (2) the demonstration at geographic scale that **changes in satellite radiance are not equivalent to changes in biologically modelled celestial polarization cues for both spiders and moths**. Avoid saying a spectral correction yields exact light emitted or perceptual detectability without validation.

*Search conducted 20 October 2026. This is an evidence synthesis for study planning, not a systematic review. Source links are stable URLs where possible; figures and years are reported as sources define them.*
