import importlib.util,json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
spec=importlib.util.spec_from_file_location('runner',Path(__file__).resolve().parents[1]/'src/register_spider_raw_batch.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class RunnerTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
  model=self.root/'model/hps';model.mkdir(parents=True);(model/'provenance.json').write_text('{}')
  self.a=SimpleNamespace(raw_batch=self.root/'batch.json',batch=self.root/'old',model_run=model.parent,output=self.root/'out',sites=None,moon_coordinates=None)
  self.rows=[]
  for site in ['museum','wells']:
   f=self.root/(site+'.json');f.write_text('{}');self.rows.append(dict(site=site,status='prepared_not_registered',returncode=0,manifest_path=str(f),manifest_sha256=m.digest(f)))
  self.save()
 def save(self):self.a.raw_batch.write_text(json.dumps(dict(status='conversion_complete_not_registered',camera_channel='green',sites=self.rows)))
 def test_hash_rejection_before_output(self):
  self.rows[1]['manifest_sha256']='bad';self.save()
  with self.assertRaisesRegex(ValueError,'manifest changed'):m.run_batch(self.a)
  self.assertFalse(self.a.output.exists())
 def test_missing_geometry_reported_with_review(self):
  s=m.run_batch(self.a);self.assertEqual(len(s['sites']),2)
  self.assertTrue(all(r['status']=='pending_inputs_or_registration' for r in s['sites']))
  page=(self.a.output/'review.html').read_text();self.assertIn('Download review JSON',page);self.assertIn('wells',page)
 def test_child_pending_and_failed_process_continue(self):
  for r in self.rows:
   f=self.a.batch/r['site']/'registration/rotation_analysis.json';f.parent.mkdir(parents=True);f.write_text('{}')
  def fake(command,**kwargs):
   site=command[command.index('--sites')+1];out=Path(command[command.index('--output')+1]);out.mkdir(parents=True)
   (out/'locked_site_references.json').write_text(json.dumps({'sites':{site:dict(status='pending_registration_not_for_comparison',transforms=[])}}))
   return SimpleNamespace(returncode=1 if site=='museum' else 0)
  with patch.object(m.subprocess,'run',side_effect=fake):s=m.run_batch(self.a)
  self.assertEqual(s['sites'][0]['status'],'registration_process_failed');self.assertEqual(s['sites'][1]['status'],'pending_registration_not_for_comparison')
if __name__=='__main__':unittest.main()
