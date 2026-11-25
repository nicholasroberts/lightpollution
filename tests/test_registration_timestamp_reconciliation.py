"""Timestamp reconciliation keeps legacy and accepted RAW metadata distinct."""
import ast
from pathlib import Path
import unittest

path=Path(__file__).resolve().parents[1]/'src/lock_spider_site_registration.py'
node=next(n for n in ast.parse(path.read_text()).body
    if isinstance(n,ast.FunctionDef) and n.name=='reconcile_raw_timestamps')
namespace={'Path':Path}
exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),namespace)
reconcile=namespace['reconcile_raw_timestamps']


class TimestampReconciliationTests(unittest.TestCase):
    def test_explicit_mapping_preserves_both_times(self):
        frames={1:dict(source_path='DSC_2259.NEF',datetime_original='2018:08:21 21:23:49')}
        old={1:'2018:08:21 22:37:47'}
        mapping=[dict(exposure=1,stem='DSC_2259',datetime_original='2018:08:21 21:23:49')]
        with self.assertRaises(ValueError):reconcile(frames,old)
        result=reconcile(frames,old,mapping)
        self.assertEqual(result[0]['old_tiff_datetime_original'],old[1])
        self.assertEqual(result[0]['raw_datetime_original'],frames[1]['datetime_original'])
        with self.assertRaises(ValueError):
            reconcile(frames,old,[dict(mapping[0],stem='wrong')])
        with self.assertRaises(ValueError):
            reconcile(frames,old,[dict(mapping[0],datetime_original='wrong')])

    def test_matching_times_need_no_exception(self):
        frame=dict(source_path='x.NEF',datetime_original='2018:08:21 21:23:49')
        self.assertEqual(reconcile({1:frame},{1:frame['datetime_original']}),[])


if __name__=='__main__':unittest.main()
