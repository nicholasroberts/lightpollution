"""Joint natural-pattern fits for the two previously reviewed Oxford sites."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import numpy as np

SITES = ('WClawn', 'oxfordbussinesspark')

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def prepare_inputs(root):
    """Check the accepted signed data lineage before writing any outputs."""
    lock_path = root/'registration/locked_site_references.json'
    locks = json.loads(lock_path.read_text())['sites']
    checks = json.loads((lock_path.parent/'reference_sha256.json').read_text())
    prepared = []
    for site in SITES:
        path = lock_path.parent/site/'registration/linear_stokes.npz'
        for source in (lock_path, path):
            key = str(source.relative_to(lock_path.parent))
            if checks.get(key) != digest(source):
                raise ValueError('Saved registration checksum differs: '+str(source))
        lock = locks[site]
        if lock.get('radiometry_status') != 'black_subtracted_linear_camera_green':
            raise ValueError('Require accepted linear camera-green Stokes: '+site)
        convention = lock.get('analyzer_stokes_convention')
        if not convention or convention.get('applied_U_multiplier') not in (-1, 1):
            raise ValueError('Require recorded signed analyzer convention: '+site)
        if lock.get('moon_pixel_xy') is None or not lock.get('astronomical_anchor'):
            raise ValueError('Require reviewed observed Moon and astronomical anchor: '+site)
        pilot_path = root/site/'fit/twilight_moon_fit.json'
        pilot = json.loads(pilot_path.read_text())
        if (pilot['source_lock_sha256'] != digest(lock_path)
                or pilot['stokes_sha256'] != digest(path)
                or pilot['camera_channel'] != 'green'
                or pilot['handedness_hypothesis'] not in ('original', 'reversed')):
            raise ValueError('Prior reviewed fit does not match signed inputs: '+site)
        with np.load(path, allow_pickle=False) as z:
            arrays = {k:np.array(z[k]) for k in ('I', 'Q', 'U')}
        if arrays['I'].ndim != 2 or any(v.shape != arrays['I'].shape for v in arrays.values()):
            raise ValueError('Stokes shapes differ: '+site)
        prepared.append((site, path, pilot_path, pilot, arrays))
    return lock_path, prepared

def export_reused_stokes(output, site, path, lock_path, arrays):
    """Adapt saved signed camera Stokes to the north-rotation interface unchanged."""
    output.mkdir(parents=True)
    np.savez_compressed(output/'measured_stokes_and_polarization.npz',
                        I_pair=arrays['I'], Q=arrays['Q'], U=arrays['U'])
    report = dict(site=site, source_lock_sha256=digest(lock_path),
                  source_stokes_sha256=digest(path), source_stokes_path=str(path),
                  status='reuse_of_checksum_verified_signed_camera_stokes',
                  raw_recomputed=False, registration_refitted=False,
                  additional_analyzer_sign_correction=False,
                  note='Interface export only; I_pair is saved I. I/Q/U copied unchanged.')
    (output/'stokes_audit.json').write_text(json.dumps(report, indent=2)+'\n')

def render(root, output, *, supplemental=False, sites=SITES, title=None, stem=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from plot_paper1_two_site_solar_lunar import render_saved_site
    plt.rcParams.update({'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig, axes = plt.subplots(4, 4, figsize=(14,15), layout='constrained')
    reports = []
    for index, site in enumerate(sites):
        base = root/site
        config = SimpleNamespace(fit=base/'fit', north=base/'north',
            reference_lock=Path((base/'source_lock_path.txt').read_text().strip()))
        report = render_saved_site(config, fig, axes[2*index:2*index+2,:])
        if report['site'] != site:
            raise ValueError('Unexpected site in saved fit')
        for row, quantity in ((2*index,'DoLP'), (2*index+1,'AoP')):
            axes[row,0].text(-.04,.5,site+' — '+quantity,
                transform=axes[row,0].transAxes, rotation=90,
                ha='right',va='center',fontsize=11)
        reports.append(report)
    fig.suptitle((title if title is not None else ('Supplementary Figure 2 — ' if supplemental else '')+
        
        'WClawn and Oxford Business Park — joint solar/lunar fits\n'
        'Shared mixture and separation within each site; no ALAN or fitted dilution'),fontsize=13)
    output.mkdir(parents=True)
    for ext in ('png','pdf','svg'):
        fig.savefig(output/((stem or ('supp_fig02' if supplemental else 'solar_lunar_two_sites'))+'.'+ext),dpi=300)
    plt.close(fig)
    (output/'figure_values.json').write_text(json.dumps(reports,indent=2,allow_nan=False)+'\n')

def run(a):
    if a.output.exists() or a.figure_copy.exists():
        raise ValueError('Use fresh output and figure-copy directories')
    lock_path, prepared = prepare_inputs(a.checked_root)
    rows = []
    for site, path, pilot_path, pilot, arrays in prepared:
        base = a.output/site
        export_reused_stokes(base/'camera_export', site, path, lock_path, arrays)
        (base/'source_lock_path.txt').write_text(str(lock_path.resolve())+'\n')
        def call(script, *args):
            subprocess.run([sys.executable,str(Path(__file__).with_name(script)),*map(str,args)],check=True)
        call('rotate_spider_measured_north.py',
             '--reference-lock',lock_path,'--audit',base/'camera_export',
             '--site',site,'--handedness',pilot['handedness_hypothesis'],
             '--output',base/'north')
        call('fit_paper1_solar_lunar_angle.py',
             '--reference-lock',lock_path,'--north',base/'north','--site',site,
             '--neutral-separation-deg',70,'--fit-neutral-separation',
             '--output',base/'fit')
        fit = json.loads((base/'fit/solar_lunar_fit.json').read_text())
        row = dict(site=site, utc=fit['solar']['utc'],
            sun_elevation_deg=fit['solar']['sun_elevation_deg'],
            sun_compass_azimuth_deg=fit['solar']['sun_compass_azimuth_deg'],
            moon_elevation_deg=fit['astronomical_anchor']['moon_elevation_deg'],
            moon_compass_azimuth_deg=fit['astronomical_anchor']['moon_compass_azimuth_deg'],
            solar_fraction=fit['solar_fraction'], lunar_fraction=fit['lunar_fraction'],
            neutral_separation_deg=fit['neutral_separation_deg'],
            fraction_on_boundary=fit['fraction_on_boundary'],
            separation_on_boundary=fit['separation_on_boundary'])
        for model in ('lunar_only','solar_only','combined'):
            row[model+'_withheld_rms_deg'] = fit['metrics'][model]['pooled_withheld']['rms_deg']
        rows.append(row)
        (base/'reuse_provenance.json').write_text(json.dumps(dict(
            prior_fit_path=str(pilot_path),prior_fit_sha256=digest(pilot_path),
            retained_handedness=pilot['handedness_hypothesis'],
            old_DoLP_mixture_reused=False,initial_separation_deg=70),indent=2)+'\n')
    (a.output/'site_fit_summary.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    with (a.output/'site_fit_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    render(a.output,a.figure_copy)
    print(json.dumps(rows,indent=2))
    print('Saved comparison:',a.figure_copy)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checked-root',type=Path,default=Path('data/processed/spider_sign_checked_01'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--figure-copy',type=Path,required=True)
    a=p.parse_args()
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))

if __name__=='__main__':main()
