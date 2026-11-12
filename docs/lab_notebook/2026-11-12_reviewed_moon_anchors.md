# Reviewed Moon anchors — 2026-11-12

Investigator reviewed ephemeris intensity overlay as 'pretty good'. Candidate source regions now coincide with visible Moon regions. Cottage original EXIF DateTimeOriginal 2018-08-21 22:37:47 Europe/Stockholm ->20:37:47 UTC; elevation11.8862deg, compass azimuth194.3726deg, candidate[1586,613]. Reference field date2018-09-21 differs by one month: preserve flag; no reference-date rewrite. Oxford original2019-05-15 00:15:50 BST ->2019-05-14 23:15:50 UTC; elevation32.4313deg, azimuth224.2925deg, candidate[580,702]. Midnight date rollover is consistent with May14 observing night; not by itself evidence of metadata error.

Added --moon-anchors to compare_wang2016_measured_skies.py. Verify anchor source-report hash. Use ephemeris elevation for source radial position under level-camera equisolid mapping; candidate image direction supplies bearing. No free source fit, global north assumption or matching of historical PRC. Review remains approximate because actual camera tilt and source centroid are not independently calibrated. Date flags retained. Constant AoP offset still fitted separately; no DoLP amplitude/ALAN amount fitted.

Existing two comparison tests pass; new direct check confirms ephemeris radius with image bearing. Real revised sky comparison remains on EcoVis.

Run prior comparison command with --moon-anchors data/processed/moon_cottage_oxford_anchors_01/moon_ephemeris_anchors.json and fresh output data/processed/wang_cottage_oxford_comparison_02. Upload comparison PNG and JSON. No new ephemeris or registration search required.
