# Figure 1A — measured Cottage and museum

The complete subfigure A has two rows (Cottage, museum) and two columns (measured green DoLP, measured green AoP). One A label applies to the entire assembly. Use north up, east left, white 60-degree full-field circle, black observed-Moon marker, linear viridis DoLP 0–1 and cyclic romaO AoP -90–90 degrees. Absolute analyzer zeros remain uncalibrated.

Read existing north-up arrays, with no new rotation, interpolation, registration or model fit. Verify north/audit/reference-lock checksums before use. Calculate reported field means and sample spatial SD from original camera I_pair,Q,U, before north resampling, using the same valid-pixel definition as the single-site panel.

Accepted recorded observations: Cottage 0.3039018694927643 +/- 0.09725326200742362, 180694 included and 165 excluded; museum 0.139762 +/- 0.040131, 180850 included and 9 excluded (museum values rounded as previously reported). Fresh script output records full precision, not hard-coded values.

## EcoVis command

```bash
python src/plot_paper1_fig01a.py \
  --output data/processed/paper1_fig01a_combined_01 \
  --figure-copy figures/paper1/fig01a/v01
```

Defaults select Cottage paper1_cottage_north_01/{north,audit/Cottage,registration/locked_site_references.json}; museum paper1_museum_north_01, paper1_museum_stokes_audit_01/museum and spider_registration_green_all_sites_01/runs/museum/locked_site_references.json. Override the corresponding CLI paths if products are relocated; source checksum matching is required.

Both directories must be new. The versioned figures copy contains PNG/PDF/SVG, figure_values.json, field_statistics.csv and archive_manifest.json with verified SHA256 of copied files. Existing versions are never overwritten. This is overwrite protection, not a filesystem security lock. Keep processed source arrays and their provenance; artwork snapshots do not replace the scientific data. Generated figure files remain excluded from Git by the existing figures/paper1/.gitignore; commit scripts and labbook, and back up generated snapshots with research data.

For other completed figure-only output folders, use:

```bash
python src/paper1_figure_archive.py \
  --source PATH_TO_COMPLETED_FIGURE_OUTPUT \
  --destination figures/paper1/FIGURE_ID/VERSION
```

The helper copies the entire supplied folder, so select the figure-output folder rather than a large RAW or reconstruction tree. Future figure runners can call archive_figure(output, destination) after calculation.

## Timing correction required before freezing 1B

The displayed old Cottage 1B uses TIFF time 2018-08-21 22:37:47 Europe/Stockholm (20:37:47 UTC). Accepted Cottage RAW reference uses 21:23:49 with explicit +01:00 offset (20:23:49 UTC), 13 min 58 s earlier. Its Moon azimuth/elevation are 191.15960520703004 / 12.31370836118817 degrees. Regenerate theory from the new lock before freezing an updated 1B; do not merely edit the old title. This combined 1A assembly does not recalculate 1B or claim the old theory is time-matched.
