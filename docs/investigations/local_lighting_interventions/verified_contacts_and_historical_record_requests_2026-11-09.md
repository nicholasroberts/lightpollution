# Stage A — verified contacts and requests for historical lighting records

**Status: prepared 9 November 2026; NOT sent.** These are targeted requests for original environmental, technical and geospatial records to define 2014–2024 interventions independently of VIIRS. The purpose is to establish treated and untreated lighting footprints within the **fixed NUTS2016** containers, not to modify any locked VIIRS/PRC products.

## Findings from additional searches

| Area | Contact (independently verified in publisher's material) | Already public | Not found publicly |
|---|---|---|---|
| Échirolles FRK24 | **courrier@echirolles.fr**, attention *Direction des services techniques / service éclairage public*; technical services tel **+33 4 76 20 64 01**. Older directory also lists courrier@ville-echirolles.fr; preferential current general contact is the shorter @echirolles.fr. | 2022 municipal notices; 2023 metropolitan slide with **EP CITE poste 079** metered reduction Oct 2022; named always-on road exceptions and LED share 5% (2022), 40% (2024). | 2022 geocoded circuits/luminaires and switching dates, actual EP079 feeder footprint, per-lamp replacement register. |
| Coventry UKG33 | **PFIStreetLighting@coventry.gov.uk** (streetlighting team); **infogov@coventry.gov.uk** (formal information request). | Official **51-page street/ward PDF and XLSX**, but those do not expose asset IDs or the exact May–Jul 2024 individual switch-off dates. | 2024 as-switched column map; exact schedule by column; street portions and exemptions; as-at-2024 lamp inventory. |
| Portsmouth UKJ31 | **foi@portsmouthcc.gov.uk** (council information); **asset.management@portsmouth.gov.uk** (general asset register enquiries, may not hold lighting operational logs). | Historical 2017 programme decision, a *zone status* map on old council PDF, 2019 pilot accounts, **4 Feb 2020** adaptive-lighting cabinet approval and council's 2022 guidance. | Dated per-column 2018–21 conversion history, pilot control-place identity, citywide CMS dimming command/commissioning dates. |
| Thanet UKJ44 | **freedomofinformation@kent.gov.uk** (**Kent County Council**, not Thanet District Council). | Current Kent lighting coordinates/asset IDs through district CSV and GIS point layer; 2016/17 and 2018/19 phase descriptions. | 2015–19 asset dated conversion and actual return to all-night lighting, including delays for concrete columns. |
| Bristol UKK11 | **lighting@bristol.gov.uk** (streetlighting); **foi@bristol.gov.uk** (formal request, if necessary). | Current council streetlight point dataset, asset IDs and CMS flag; 2023 start and subsequent project summaries; vendor description of 38 schedules and half-hourly telemetry. | Historical 2023–24 dates of lantern swaps versus CMS-only activation, pre-2023 lamp types and actual 2023–24 operational profiles. |

**Caution Échirolles:** December 2023 PDF labels some action dates as *2023* in its table while the same slide's energy graph unmistakably identifies the October **2022** drop and official council announcements made in 2022 independently support that date. Request the original actuation/change log and clarification. The graph is energy use, not luminous flux.

**Caution Portsmouth:** the zone map establishes a historical rollout *status snapshot*, not verified per-lamp installation dates. The unmodified pilot-control neighbourhood reported in contemporary press remains unidentified.

**Caution Coventry:** the public street list lacks per-lamp coordinates, days and potentially includes only portions of some streets. Do not automatically treat all lamp columns on each listed road.

## Suggested exact requests — ready to send

### 1. Échirolles — French original

**To:** courrier@echirolles.fr  
**Attention:** Direction des services techniques / service éclairage public  
**Subject:** Demande de données historiques — extinction de l'éclairage public (2022) et rénovation LED (2022–2024)

Madame, Monsieur,

Je mène à l'Université de Bristol une recherche scientifique sur l'évolution de l'éclairage nocturne entre 2014 et 2024, à partir des observations satellitaires VIIRS. Le programme d'extinction engagé à Échirolles en 2022 constitue un cas d'étude particulièrement intéressant. Pour définir les secteurs concernés indépendamment des données satellitaires, pourriez-vous transmettre, dans la mesure des documents disponibles :

1. L'inventaire géolocalisé des points lumineux communaux en **2022** (identifiant, coordonnées ou adresse/voie, type de lampe, puissance et, si disponible, flux lumineux), ainsi que les principales modifications **2022–2024**.
2. Un plan ou une liste des **voies, circuits ou armoires d'éclairage effectivement éteints**, avec les dates de mise en application et horaires à partir d'**août et octobre 2022**, et la liste des voies restées éclairées (notamment les grands axes).
3. Pour **« EP CITE poste 079 »**, présenté lors du forum métropolitain de décembre 2023 : le secteur ou les luminaires alimentés et, si disponibles, les **consommations électriques mensuelles ou journalières 2020–2023** ; un agrégat par secteur convient.
4. Les dates et secteurs de remplacement des anciennes lampes par des **LED en 2022–2024**, avec la technologie d'origine, la puissance, les caractéristiques des nouveaux luminaires et les changements d'horaires éventuels.
5. Une clarification de la chronologie : le graphique du document métropolitain situe la rupture de consommation en **octobre 2022**, alors que son tableau indique **2023** pour certaines mesures.

Un tableau CSV/Excel et des données géographiques GeoJSON/Shapefile seraient idéaux, mais des plans, listes de rues ou copies des documents existants seraient déjà très utiles. Il ne nous faut aucune donnée personnelle ni schéma électrique sensible. Pourriez-vous transmettre ce message au service compétent si nécessaire ?

Je vous remercie pour votre aide et serai heureux de citer les sources et contributions de la Ville.

Bien cordialement,

[Nom]  
School of Biological Sciences, University of Bristol  
[Adresse électronique]

### 2. Coventry — English

**To:** PFIStreetLighting@coventry.gov.uk  
**Cc (optional for formal request):** infogov@coventry.gov.uk  
**Subject:** Research data request: actual 2024 part-night streetlight commissioning and asset locations

Dear Street Lighting Team,

I am studying historical changes in night-time lighting across Europe using VIIRS satellite observations. Coventry's May–July **2024** part-night switch-off offers a valuable independently documented intervention, before the main **2025** LED programme.

I have the council's publicly available street/ward switch-off PDF and Excel list. Would you please supply the underlying **as-implemented 2024** records, preferably as an extract from the council's/PFI contractor's asset register:

1. The **asset/column ID, point location** (British National Grid or lat/long), street and ward for lights affected by the 2024 part-night programme, identifying partially treated streets separately.
2. Each column's **actual switch-over/programming date** in May–July 2024, programmed **weekday/weekend switch-off times**, and any exceptions, reversals or operating changes to 31 December 2024.
3. Comparable **unaffected/exempt columns**, with the reason for exemption, so that appropriate within-city controls can be selected independently of satellite measurements.
4. The **lamp technology and rated wattage/output as at early 2024**, including any LED installations predating the principal programme launched in 2025, and available feeder/aggregate electricity data for 2023–24.

CSV/Excel and GeoJSON/GIS are preferred; any existing 2024 programme extract, mapping or contractor completion certificate would be helpful. No personal information is needed. Please pass the request to Balfour Beatty/the relevant PFI records holder where appropriate.

Kind regards,  
[Name]  
School of Biological Sciences, University of Bristol

### 3. Portsmouth — English

**To:** foi@portsmouthcc.gov.uk  
**Subject:** Environmental information request: Portsmouth 2018–2021 streetlight LED installation and adaptive dimming records

Dear Information Governance Team,

I am undertaking academic research into how LED replacement and subsequent lighting dimming separately affect satellite-observed night-time radiance. Could you provide recorded information relating to the Portsmouth LED rollout (**2018–2021**) and the adaptive-lighting pilot (**2019**) / citywide scheme (**summer 2020**)?

1. A **dated street/column-level LED conversion register** with asset ID, position, original and replacement lamp type, rated watts/lumens and, where available, LED colour temperature.
2. The **2019 pilot-area maps, control and treatment assignments** for Baffins, Cosham, St Jude, Milton and Fratton, including identification of the neighbourhood deliberately left unchanged and the precise trial dates.
3. **CMS dimming schedules and actual activation/commissioning dates** (and recorded output set points) for the pilot and 2020 citywide rollout, including exceptions such as main roads, city centre and heritage lighting.
4. Relevant 2018–21 **metered energy or CMS operational summaries**, ideally spatially tied to street sections or lighting groups.

CSV/Excel and GIS are ideal, but archived zone maps, decision appendices or programme extracts are useful. Please treat this as a request for recorded environmental information under the Environmental Information Regulations, or direct it to the lighting team/contractor holding the underlying data. We do not require any personal information.

Kind regards,  
[Name]  
School of Biological Sciences, University of Bristol

### 4. Thanet — English

**To:** freedomofinformation@kent.gov.uk  
**Subject:** Environmental information request: Thanet 2016–2019 LED conversion and restoration of all-night operation

Dear Kent County Council Information Governance Team,

I am comparing historical streetlight interventions with satellite-observed night-time radiance. I have located the council's current streetlight asset locations and the reported **6,922** Thanet conversions during December 2016–May 2017, with **3,009** outstanding as of 2017.

Could you supply the historical asset/programme records needed to establish when and where the changes actually occurred?

1. **2015–2019 lighting asset IDs and coordinates**, ideally using the same persistent IDs as the current council streetlight location CSV/GIS.
2. For each affected asset, the **actual date of LED conversion/commissioning**, previous and new luminaire technology and watts/lumens, with a flag for the **second-phase concrete-column and heritage delays**.
3. **Dates when night-time lighting was restored to all-night operation** at each asset or lighting group, together with the pre-conversion part-night switch-off hours and post-conversion operating times.
4. Any **monthly phase-completion summaries, contractor records, maps or 2015–19 lighting electricity data** held for Thanet.

CSV/Excel and an asset-linked GIS extract are preferred; archived reports would also be valuable. These records are requested to define intervention footprints before examining satellite outcomes, and no personal details are needed.

Kind regards,  
[Name]  
School of Biological Sciences, University of Bristol

### 5. Bristol — English

**To:** lighting@bristol.gov.uk  
**Cc (optional for formal request):** foi@bristol.gov.uk  
**Subject:** University of Bristol research request: historic 2023–2024 streetlight LED/CMS installation and operation data

Dear Bristol Street Lighting Team,

I am undertaking a University of Bristol research project comparing night-time satellite observations (VIIRS) against independently documented changes to street lighting. Bristol's rolling LED and central management system (CMS) programme, beginning in **March 2023**, would be a particularly valuable local case.

The council already publishes a streetlight GIS layer including unit identifiers and a current CMS-control flag. To distinguish **actual 2023–2024** changes from current asset status, could you provide:

1. An **asset-ID-linked installation/commissioning history** from January 2023 to December 2024: actual lantern replacement date, CMS-node-only installation date, CMS activation date and streets/wards covered.
2. The **previous and new lamp technology**, watts, nominal lumen output and CCT (or product specifications), including the important distinction between replacement LED lanterns and pre-existing LEDs that received CMS nodes only.
3. The **historical schedule profiles** and group/asset assignments effective in 2023 and 2024, including midnight/early-morning dimming times and output set points; not just the currently active 38 schedules.
4. Where readily available, **half-hourly energy/activity summaries or commissioning/exception logs** around the VIIRS overnight acquisition period, or anonymised/aggregated street-group summaries.
5. Any dates of **not-yet-converted** streets as at 31 December 2024, enabling a matched comparison within Bristol's original geographical boundary.

A CSV/Excel extract keyed to the public streetlight GIS IDs, plus any available GIS coordinates or metadata, would be ideal. No personal data is needed; if Centregreat or the CMS provider holds the records, I would appreciate a referral or an introduction. Historical data through **2024** is the priority, not the later completed scheme.

Many thanks,  
[Name]  
School of Biological Sciences, University of Bristol

## Verified links and access notes

- Échirolles municipal public lighting jurisdiction / technical services: https://www.echirolles.fr/vie-quotidienne/cadre-de-vie/voiries-espaces-publics-proprete
- Échirolles **current** general correspondence mailbox in 2024 council press release footer: https://www.echirolles.fr/sites/default/files/2024-10/20241001_CP_RepriseCasino_Uexpress%282%29.pdf
- Échirolles 2022 notice: https://www.echirolles.fr/actualites/un-nouveau-plan-de-sobriete-energetique-pour-la-ville ; metropolitan meter evidence (slide 7 in PDF, with internal date inconsistency): https://planclimat.grenoblealpesmetropole.fr/cms_viewFile.php?idtf=1784&path=Parole-aux-partenaires.pdf
- Coventry street list: https://www.coventry.gov.uk/downloads/download/7837/part-night-switch-off-alphabetical-list-of-streets-and-ward ; contact: https://www.coventry.gov.uk/street-scene/coventry-led-streetlighting-upgrade/3 ; information team: https://www.coventry.gov.uk/freedom-information-environmental-information-regulations/foi-eir-disclosure-log
- Portsmouth 2018/19 rollout-zone map: https://www.portsmouth.gov.uk/wp-content/uploads/2020/05/trv-led-streetlight-map-zone-2.pdf ; adaptive-lighting approval: https://portsmouth.moderngov.co.uk/documents/s26318/Record%20of%20decisions%20Cabinet%204%20Feb%202020.pdf ; FOI: https://www.portsmouth.gov.uk/services/council-and-democracy/transparency/freedom-of-information/publication-scheme ; asset-management general: https://www.portsmouth.gov.uk/services/council-and-democracy/transparency/freedom-of-information/lists-and-registers
- Kent streetlight points/codes (district CSV): https://www.kent.gov.uk/roads-and-travel/what-we-look-after/street-lights/faulty-street-lights ; GIS layer: https://webapps.kent.gov.uk/server/rest/services/Test/KentViewPortalTest/MapServer/44 ; FOI: https://www.kent.gov.uk/about-the-council/information-and-data/access-to-information/freedom-of-information-requests
- Bristol lighting team: https://services.bristol.gov.uk/files/documents/4350-03-8-2-street-lighting ; FOI: https://www.bristol.gov.uk/council/data-protection-and-foi/freedom-of-information-foi ; asset GIS: https://opendata.bristol.gov.uk/datasets/bcc::streetlights-3

## Request strategy and caveats

- Make a **collaborative research request to the streetlighting team first**, where direct technical contacts exist; formally ask for recorded environmental data under **EIR 2004** if substantive records cannot be supplied. UK authorities typically target a response within 20 working days, subject to exceptions.
- Each request is **time-limited and dataset-specific**; ask for existing extracts/exports and documents instead of asking staff to perform a custom analysis.
- If column-level histories cannot be shared, accept **dated street/feeder polygons, phase-completion records and aggregate electricity per feeder** as useful fallbacks. Allow exclusion of infrastructure/security-sensitive detail.
- Distinguish **actual installation, programme approval, commissioning and lighting operation**. The relevant satellite test requires the state near the VIIRS night-time overpass, not simply nominal energy consumption.
- Do not republish third-party data to public GitHub unless licensing/redistribution rights are confirmed; metadata and source URLs may be recorded.
- **No request has been sent**. Contact addresses have been confirmed on source material, but successful delivery has not been tested.
