# ChurchHill star picking and zoology acceptance — 2026-11-27

The investigator accepts zoology following inspection of the reference image
and registration readiness report. Record the accepted source precisely:
- Lock SHA256: c81554432b3fc143f270afab9ab6b130aba58c90cb9ca53e8dd375834a7d3839
- Stokes SHA256: 74f4825d60431e397041acec32a61825c3c5b4beb36ef375082ebc2729c77c4f
- Reference image SHA256: a2beb0515f5596626a2d9835c2120fd5d0be40c5adc6e35df87f186060f31879
- Accepted Moon pixel: [1420, 631].
- Nominal analyzer sequence is positive; preserve camera U.
- Thin cloud is visible near the Moon; retain this observation in interpretation.
This entry records investigator acceptance without mutating the original locks.

ChurchHill has no identified Moon. Its saved ephemeris gives elevation 8.0675
degrees and compass azimuth 126.3410 degrees, at 2019-04-20T21:47:51 UTC
(local 22:47:51, existing timezone assumption). Next identify stars independently
of polarization, using the same unchanged picker used for PennCC:

```bash
python src/pick_spider_reference_stars.py \
  --reference-lock data/processed/spider_registration_green_all_sites_01/runs/ChurchHill/locked_site_references.json \
  --site ChurchHill \
  --output data/processed/paper1_ChurchHill_star_picker_01
xdg-open data/processed/paper1_ChurchHill_star_picker_01/pick_stars.html
```

Pick the clearest stars spread across the sky, ideally four or more, then
download ChurchHill_reference_star_points.json and the marked PNG.
Labels P1, P2 etc suffice; identities will be tested against the capture
geometry afterwards. Brightness and zoom controls affect display only.
The helper preserves original pixel coordinates and source/image checksums.
No code adaptation is needed: the existing picker already supports any site.

ChurchHill's numerical registration and negative nominal analyzer sequence
pass the saved checks; the U correction is not yet applied. Do not infer
orientation from the polarization pattern. The picked stars are candidates,
not a celestial solution or a confirmed Moon position.
