"""Accepted exposure-time tolerance does not change RAW values or other settings."""
import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('raw_tolerance',Path(__file__).resolve().parents[1]/'src/prepare_spider_linear_raw.py')
raw=importlib.util.module_from_spec(spec);spec.loader.exec_module(raw)


class ExposureToleranceTests(unittest.TestCase):
    def test_small_accepted_difference_and_strict_default(self):
        ref=dict(iso=3200,exposure_seconds=59.9,f_number=3.5)
        frame=dict(ref,exposure_seconds=59.7)
        self.assertTrue(raw.exposure_settings_changes(frame,ref))
        self.assertEqual(raw.exposure_settings_changes(frame,ref,.2),[])
        self.assertEqual(frame['exposure_seconds'],59.7)
        self.assertEqual(ref['exposure_seconds'],59.9)

    def test_other_settings_and_larger_difference_still_rejected(self):
        ref=dict(iso=3200,exposure_seconds=59.9,f_number=3.5)
        for frame in [dict(ref,iso=1600),dict(ref,f_number=4),dict(ref,exposure_seconds=59.6)]:
            self.assertTrue(raw.exposure_settings_changes(frame,ref,.2))


if __name__=='__main__':unittest.main()
