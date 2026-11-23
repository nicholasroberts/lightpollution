import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src/paper1"))
import _workflow
from audit_sites import audit


class PaperWorkflowTests(unittest.TestCase):
    def test_inventory_keeps_pending_and_discovers_sites(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = root / "sites.csv"
            cfg.write_text("site,site_id,location,setup_state\nMissing,s1,Test,pending\n")
            raw = root / "raw"
            (raw / "New" / "nested").mkdir(parents=True)
            (raw / "New" / "nested" / "DSC_001.NEF").write_bytes(b"fixture")
            records = {r["site"]: r for r in audit(cfg, raw)}
            self.assertEqual(set(records), {"Missing", "New"})
            self.assertEqual(records["Missing"]["raw_directory"], "absent")
            self.assertEqual(records["New"]["nef_count"], 1)
            self.assertFalse(records["New"]["moon_anchor_recorded"])
            self.assertEqual(records["New"]["scientific_approval"], "not_assigned_by_inventory")

    def test_dispatch_dry_run_cannot_execute_forwarded_command(self):
        result = subprocess.run([sys.executable, str(ROOT / "src/paper1/run_all.py"),
                                 "--stage", "conversion", "--", "--output", "do-not-create"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("Dry-run only", result.stdout)

    def test_pending_analysis_fails_without_result(self):
        result = subprocess.run([sys.executable, str(ROOT / "src/paper1/compare_equal_luminance.py")],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("PENDING", result.stderr)

    def test_adapter_records_completed_command_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "src").mkdir()
            (root / "config").mkdir()
            (root / "config" / "model.yaml").write_text("status: template\n")
            (root / "src" / "fixture.py").write_text(
                "import pathlib,sys\npathlib.Path(sys.argv[2]).mkdir()\n")
            output = root / "result"
            with patch.object(_workflow, "ROOT", root), patch.object(_workflow, "CONFIG", root / "config"):
                with self.assertRaises(SystemExit) as ended:
                    _workflow.forward("fixture.py", ["--output", str(output)])
                self.assertEqual(ended.exception.code, 0)
                record = json.loads((root / "result.paper1-provenance.json").read_text())
                self.assertEqual(record["status"], "completed")
                self.assertEqual(record["command"][-1], str(output))
                self.assertTrue(record["source_script_sha256"])
                with self.assertRaises(SystemExit):
                    _workflow.forward("fixture.py", ["--output", str(output)])


if __name__ == "__main__":
    unittest.main()
