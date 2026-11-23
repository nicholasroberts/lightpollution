# Paper 1 workflow checks

Run from the repository root:

```bash
python -m unittest discover -s tests/paper1 -v
```

These checks exercise pending-site retention, RAW filename inventory without
pixel reads, dry-run dispatch, explicit failure of unfinished analyses and
provenance/overwrite behavior. They are not validation of the sky physics or
experimental results. Existing scientific tests remain alongside legacy code.
