# Cottage corrected full-field results — 2026-11-12

Reviewed uploaded polarization_summary(2).json and matched plots. Confirmed angular radius 30 degrees, full field 60 degrees, pixel radius 240.0361762, centre [1001.7451430,659.5240893]. All 181020 selected pixels physical in both reconstructions; zero excluded.

Intensity-retaining fixed-eye PRC: linear Stokes .4540884414; legacy DoP/AoP .4468182387. Historical reference .439: absolute differences +.0150884414 (+3.437%) and +.0078182387 (+1.781%). Previous wrong 60-degree-radius values were .3148367066 and .2927137726. Correcting FOV therefore explains most of prior PRC deficit. It does not establish exact reproduction or validate other sites.

Mean local DoLP .3135794838 direct and .3121135578 legacy (~.47% relative difference). Effective integrated intensity-weighted DoLP .2996605185 and .2950200427. Normalized-direction PRC alternatives .4704230967 and .4636175214 remain separate; do not select method solely for proximity to .439. Unweighted orientation coherence .9887735925/.9795517920 indicates much less spatial angular cancellation than previous large field.

No numerical tuning performed to obtain reference. Remaining small differences could involve historical pixel geometry, registration, replicate averaging and intensity handling; causes unresolved. Original whole-sky AoP figure provenance remains unresolved. Current strong local angular agreement in corrected FOV means that whole-sky discrepancy need not imply a large central-field contrast discrepancy.

Next: propagate corrected full-field convention into core reconstruction/configuration/documentation and repeat on another original calibration site. Preserve separate processing methods until historical details resolved. Locked reference values and fits unchanged.
