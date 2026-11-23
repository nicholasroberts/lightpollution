"""Inventory all named sites without fitting or approving scientific checks."""
import argparse
import csv
import json
import sys
from pathlib import Path
from _workflow import ROOT, provenance


def audit(config, raw_root=None, reference_lock=None, results_root=None):
    with Path(config).open(newline="") as f:
        configured = {r["site"]: r for r in csv.DictReader(f)}
    locks = json.loads(Path(reference_lock).read_text()).get("sites", {}) if reference_lock else {}
    discovered = {p.name for p in Path(raw_root).iterdir()
                  if p.is_dir() and not p.name.startswith(".")} if raw_root else set()
    result_sites = {p.parent.name for p in Path(results_root).glob("*/prc_lamp_equivalence.json")} if results_root else set()
    names = sorted(set(configured) | set(locks) | discovered | result_sites)
    records = []
    for site in names:
        prior = configured.get(site, {})
        folder = Path(raw_root) / site if raw_root else None
        raw_count = sum(1 for p in folder.rglob("*") if p.is_file() and p.suffix.lower() == ".nef") if folder and folder.is_dir() else None
        lock = locks.get(site, {})
        record = dict(site=site, site_id=prior.get("site_id", ""), location=prior.get("location", ""),
                      setup_state=prior.get("setup_state", "new_site_requires_review"),
                      raw_directory="present" if folder and folder.is_dir() else ("absent" if raw_root else "not_inspected"),
                      nef_count=raw_count,
                      registration_status=lock.get("status", "not_in_supplied_lock"),
                      moon_anchor_recorded=bool(lock.get("moon_pixel_xy") is not None and lock.get("astronomical_anchor")),
                      analyzer_convention_recorded=bool(lock.get("analyzer_stokes_convention")),
                      equivalence_report="present_not_revalidated" if site in result_sites else "not_found",
                      scientific_approval="not_assigned_by_inventory",
                      next_action="reconcile source checks and prior visual reviews")
        records.append(record)
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", type=Path, default=ROOT / "config/paper1/sites.csv")
    p.add_argument("--raw-root", type=Path)
    p.add_argument("--reference-lock", type=Path)
    p.add_argument("--results-root", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        p.error("Use a fresh output directory")
    for path in [a.sites, a.raw_root, a.reference_lock, a.results_root]:
        if path is not None and not path.exists():
            p.error("Missing input: " + str(path))
    records = audit(a.sites, a.raw_root, a.reference_lock, a.results_root)
    if not records:
        p.error("Empty site inventory")
    meta = provenance(Path(__file__).resolve(), sys.argv[1:])
    meta["input_paths"] = {k: str(v) for k,v in vars(a).items() if v is not None}
    from _workflow import digest
    meta["input_hashes"] = {str(p): digest(p) for p in [a.sites, a.reference_lock] if p is not None}
    meta["raw_images_opened"] = False
    meta["scientific_checks_approved"] = False
    a.output.mkdir(parents=True)
    with (a.output / "site_status.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (a.output / "site_status.json").write_text(json.dumps(records, indent=2) + "\n")
    (a.output / "provenance.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Inventoried {len(records)} sites; no registration or Moon checks approved.")


if __name__ == "__main__":
    main()
