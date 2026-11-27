# ChurchHill four-star candidate solution — 2026-11-27

Read the attached ChurchHill_reference_star_points.json. Use fixed saved
reference centre/radius and the existing capture UTC 2019-04-20T21:47:51,
latitude 52.555609, longitude -2.156873, pressure=0, height=0.
A PyEphem 4.2.1 named-star search tested distinct above-horizon stars within
5 degrees of each image elevation, both bearing parities and one circular-mean
azimuth offset. No polarization objective or image rotation fit was used.

| Point | Candidate star | Position residual |
|---|---|---:|
| P1 | Arcturus | 1.664° |
| P2 | Pollux | 0.847° |
| P3 | Regulus | 0.244° |
| P4 | Procyon | 0.697° |

Best positional RMS: 1.004 degrees; reversed handedness, azimuth =
-pixel_angle -36.9111 degrees (modulo 360), x right/y down.
The next named-star assignment substitutes Castor for P2 and has RMS 2.530
degrees; candidate identification is not yet investigator accepted.
No independent holdout star, projection/tilt fit or timestamp verification.

The solution predicts Moon pixel [419.7163, 483.8928], elevation 8.0674 degrees,
azimuth 126.3420 degrees. This falls at the left-side building/skyline region
in the supplied reference image and is consistent with an obscured Moon.
It is a prediction from stars and ephemeris, not an observed Moon detection.
Do not yet mutate the accepted registration or fit a north-up pattern until
the star solution is reviewed. Points, source hashes and candidate solution
are saved in config/paper1/ChurchHill_reference_star_points.json.
