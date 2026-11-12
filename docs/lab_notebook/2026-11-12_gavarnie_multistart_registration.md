# Gavarnie registration retry — 2026-11-12

Retrieved investigator-uploaded rotation_analysis(2).json and overlay record. Report identifies image 12 at 174.065 degrees for nominal 135, inner structural correlation .0531, translation-limit flag. Image 8 initial polar estimate 73.652, refined 87.651, inner correlation .1568 and rotation-limit flag. Images 10/11 also hit angular search limits; images 7/9 recover approximately 90 degrees. Prior direct PRC .246717 and legacy .254288 versus reference .470 cannot establish reproduction under these alignment failures.

Added opt-in --multistart to analyse_spider_sky_rotation.py. For each exposure, compare refinement from polar estimate and both signed nominal camera-angle seeds, dropping near-duplicate seeds. Select by inner-sky structural correlation, never contrast or proximity to historical PRC. Save trial seeds, rolls and scores in report. Eye FOV remains 30deg radius; geometry and .4/.04 absorbances unchanged. Default mode and historical batch outputs preserved.

New synthetic regression: imposed 135deg scene roll and [-4,+5]px source displacement, intentionally incorrect initial polar roll 174deg. Retry recovers 134.9997deg, [-4.0005,+5.0002]px with correlation .9849; wrong seed remains near .066. Regression executed successfully. Does not claim actual Gavarnie correction before EcoVis rerun.

Run original batch Gavarnie geometry as input, original TIFF folder, site-id s12, --multistart, and fresh output data/processed/spider_prc_gavarnie_alignment_02. Review rotation_overlay.png, rotation_analysis.json and both method preview result.json files. No need to rerun all sites until this targeted fix is assessed. Raw images and locked reference PRC unchanged.
