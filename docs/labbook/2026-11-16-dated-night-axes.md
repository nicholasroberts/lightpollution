# Dated local-time axes and Cottage's truncated lunar interval

Cottage capture is21Aug2018 at22:37:47CEST. Rechecking the saved ephemeris
configuration at10min steps places Moon elevation +0.14deg at01:20CEST on
22Aug and -0.92deg at01:30. The model stops when Moon elevation falls to
0.1deg, approximately01:20–01:21 local. The Sun remains below-6deg until
roughly05:08CEST, so the remaining grey interval is chiefly Moon-below-horizon,
not dawn. The first valid dark interval begins approximately21:05 on21Aug.
Times here are approximate, based on the builtin ephemeris, no refraction,
nominal site coordinates and sampled crossings.

The model deliberately calculates the lunar cue only for Sun elevation<-6deg
and Moon elevation>0.1deg and<89.9deg, with further near-Moon field exclusions.
It does not model moonless starlight/airglow polarization or assign a zero PRC
to excluded times. Grey therefore means excluded/unmodelled, not no darkness
or an observed zero signal.

Nightly heatmaps, PRC heading-range plots and spectral comparison time-series
now use actual Matplotlib datetime coordinates. Tick labels show day/month/year
and24h clock time in Europe/Stockholm for Cottage and Europe/London for Oxford.
UTC calculations and CSV timestamps are unchanged. Midnight shows the next
calendar date, replacing the24,25,...hour notation. Panel titles are split over
two lines, and the spectrum label is included in the nightly figure heading.

Changes apply when figures are next generated with the standard runner; old
PNG files are unchanged. Two tests check local midnight rollover and rejection
of timezone-naive inputs. A coarse overnight plotting run checks the dated
heatmaps and ranges without changing the physical calculation.
