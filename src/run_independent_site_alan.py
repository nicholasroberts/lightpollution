"""Run independent site-ALAN sensitivity scenarios and measured-map comparisons.
No DoLP/PRC fitting. Upward-radiance contrast is only a provisional skyglow proxy.
"""
from pathlib import Path
import argparse,subprocess,sys,json
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['batch','spectra','moon-anchors','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--levels',type=float,nargs='+',default=[.3,1.,3.,10.])
    p.add_argument('--atmospheric-factor',type=float,default=.8)
    p.add_argument('--skyglow-ratio',type=float,help='Override proxy only with independently established Oxford/Cottage ground-sky ratio')
    p.add_argument('--rsr-file',type=Path)
    p.add_argument('--cottage-target-dop',type=float,help='Explicit user scenario adjustment; mean HPS field DoLP at target level')
    p.add_argument('--cottage-target-level',type=float,default=5.)
    a=p.parse_args()
    import numpy as np
    if not a.levels or any(not np.isfinite(v) or v<0 for v in a.levels) or len(set(a.levels))!=len(a.levels):p.error('Distinct finite nonnegative levels required')
    if not np.isfinite(a.atmospheric_factor) or not 0<=a.atmospheric_factor<=1:p.error('Atmosphere must be0..1')
    if a.output.exists():p.error('Use fresh output directory')
    if not a.batch.is_dir():p.error('Saved measured batch directory required')
    a.output.mkdir(parents=True);out=a.output.resolve()
    def run(script,*args):
        command=[sys.executable,str(ROOT/'src'/script),*map(str,args)]
        print('Running '+script,flush=True);subprocess.run(command,check=True)
    prep=['--spectra',a.spectra.resolve(),'--output',out/'alan_configs']
    if a.skyglow_ratio is not None:prep+=['--skyglow-ratio',a.skyglow_ratio]
    if a.rsr_file:prep+=['--rsr-file',a.rsr_file.resolve()]
    run('prepare_site_alan_scenarios.py',*prep)
    if a.cottage_target_dop is not None:
        run('adjust_cottage_alan.py','--configs',out/'alan_configs','--spectra',a.spectra.resolve(),'--moon-anchors',a.moon_anchors.resolve(),
            '--target-dop',a.cottage_target_dop,'--target-level',a.cottage_target_level,'--atmospheric-factor',a.atmospheric_factor)
    for source,label in [('HPS','hps'),('LED_mid','led')]:
        run('model_spider_nightly_prc.py','--moon-anchors',a.moon_anchors.resolve(),'--spectra',a.spectra.resolve()/(source+'_visible_input.csv'),
            '--site-alan-config',out/'alan_configs'/(source+'_site_alan.json'),'--atmospheric-factors',a.atmospheric_factor,'--alan-ratios',*a.levels,'--tilts',0,'--output',out/label)
    for level in a.levels:
        run('plot_spider_spectral_pilot.py','--spectra',a.spectra.resolve(),'--hps',out/'hps','--led',out/'led','--batch',a.batch.resolve(),
            '--atmospheric-factor',a.atmospheric_factor,'--alan-ratio',level,'--output',out/f'comparison_level_{level:g}')
    (out/'scenario_notes.json').write_text(json.dumps(dict(status='USER_ADJUSTED_COTTAGE_SCENARIO' if a.cottage_target_dop is not None else 'NO_DOP_FIT',cottage_target_dop=a.cottage_target_dop,cottage_target_level=a.cottage_target_level,levels=a.levels,meaning='Oxford HPS ALAN450–750nm / Cottage capture zenith natural intensity; Cottage follows explicit target when requested, otherwise independent site ratio; LED scaled for equal DNB-overlap source response',caveat='When cottage_target_dop is set the prior site-ratio proxy is superseded for Cottage; Oxford remains unchanged. Absolute ALAN levels are sensitivity scenarios. Upward radiance is not a measured downward skyglow ratio.'),indent=2)+'\n')
    print(f'Finished. Comparison figures are in {out}/comparison_level_*/',flush=True)

if __name__=='__main__':main()
