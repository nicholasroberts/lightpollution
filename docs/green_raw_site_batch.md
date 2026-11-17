# Audited green RAW site batch

The user confirmed retained filename/time order as the intended analyzer acquisition sequence, irrespective of gaps caused by deleting transient-light frames. Twelve retained files map to four angle groups of three repeats and quartets [1,4,7,10], [2,5,8,11], [3,6,9,12]. Waitrose's four retained files map to one [1,2,3,4] quartet. Mapping indices are retained even when conversion exclusions are applied.

`prepare_spider_raw_batch.py` reads the metadata audit, resolves folder aliases ChurcHill -> ChurchHill and WC_pitch -> WCpitch, and preserves wells as Wells Cathedral, a separate site. It validates all audited source sizes, timestamps and exposure settings against the local NEFs before decoding. Inventory changes and unreviewed settings changes stop preflight. No source file is modified. SHA256 input hashes are then calculated by the existing converter.

Two explicitly recorded exclusions avoid applying exposure normalization: museum DSC_5251 (sequence position 7, 9.4 versus 9.9s), and slabhouse DSC_5404 (position 10, 24.5 versus 24.9s). Remaining complete quartets are retained; other frames remain available as geometric references. The attached 172-file audit yields 15 sites, 170 converted TIFFs and 41 candidate complete quartets. These are acquisition groups, not validated polarimetry results.

```bash
cd /run/media/bznwr/data/lightpollution
git pull origin publication-figures

python src/prepare_spider_raw_batch.py \
  --audit data/processed/nef_metadata_audit.json \
  --raw-root data/raw/skyimages \
  --channel green \
  --output data/processed/spider_linear_green_all_sites_01
```

The output contains `mappings/SITE.json`, each site's linear TIFFs and `linear_raw_manifest.json`, individual conversion logs, and `batch_summary.json`. Conversions are sequential, one site at a time. Failed sites are recorded and the batch continues to remaining sites; the command exits nonzero if any fail. Existing output directories are never overwritten. For selected sites, add `--sites WClawn WCpitch` using canonical names. `--plan-only` validates local files and writes mapping/plan records without decoding; use a separate fresh output directory for that mode.

Upload `batch_summary.json` after running. Check failure logs where indicated. Conversion does not perform registration, Moon localization, astronomical calculations or model fitting. The previously completed museum results remain in their original folders. The four-frame converter support is preparation only; downstream multi-site registration must explicitly support that grouping before Waitrose reconstruction. Wells coordinates and reference/batch identity must be supplied independently before its astronomical model is attempted.
