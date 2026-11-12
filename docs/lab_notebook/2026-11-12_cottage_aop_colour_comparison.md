# Original-legend AoP colour comparison — 2026-11-12

User requested matching the original rainbow AoP colour scale. Recovered R source palette differs from uploaded PNG; therefore use an explicitly approximate colour progression sampled from the uploaded original PNG legend, not a claim to recover its exact historical palette. Low-to-high red/orange/yellow/green/cyan/blue/violet/magenta/red-pink. Shared explicit limits -90 to +90 degrees on both reconstructed methods; original PNG numeric range/units remain unverified. No angle rescaling, offsets, sign changes, wrapping changes, Stokes or PRC changes.

Updated audit_spider_polarization.py to apply this palette to field panels and create matched_aop.png with both complete reconstructed AoP fields in transposed original R axis order, common degree-labelled scale, equal pixel aspect and 60-degree contour. This permits comparison with full-sky original rather than only a central masked crop. Original plotted aspect ratio is not evidence for source-pixel anisotropy. Matched-polarization display retained. Syntax compilation passed; existing scientific calculations unchanged.

Run audit again using fresh output data/processed/spider_prc_cottage_polarization_audit_02. No registration rerun needed. Review matched_aop.png; identical colours alone do not validate AoP numerically.
