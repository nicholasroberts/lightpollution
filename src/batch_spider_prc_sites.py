"""Process the 18 historical spider sites; preserve per-site diagnostics and logs.

Nominal geometry is transferred from Cottage with scale based on short image
side; it is an explicit approximation, not an independent calibration per site.
Different aspect ratios/crops and failed circle fits are flagged for review.
Cottage reuses its reviewed per-exposure geometry. No additional folders enter
historical comparisons (Bristol, Sandhammeren and Wells are excluded).
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import subprocess
import sys
import numpy as np
import tifffile
from reconstruct_spider_eye_contrast import digest

ROOT=Path(__file__).resolve().parents[1]
SITES={
 'ChurchHill':('s07',['ChurchHill','ChurcHill']), 'Cottage':('s14',['Cottage']),
 'EastHorrington':('s13',['EastHorrington']), 'Gavarnie':('s12',['Gavarnie']),
 'minifactory':('s05',['minifactory']), 'museum':('s15',['museum']),
 'opensquare':('s16',['opensquare']), 'Ounsdale':('s09',['Ounsdale']),
 'oxfordbussinesspark':('s02',['oxfordbussinesspark']), 'PennCC':('s08',['PennCC']),
 'Skillinge':('s10',['Skillinge','skillinge']), 'slabhouse':('s17',['slabhouse']),
 'southhorrington':('s11',['southhorrington']), 'stannes':('s04',['stannes']),
 'Waitrose':('s03',['Waitrose']), 'WClawn':('s18',['WClawn']),
 'WCpitch':('s01',['WCpitch','WC_pitch']), 'zoology':('s06',['zoology'])}


def locate(root,aliases):
    candidates=[p for p in root.iterdir() if p.is_dir() and p.name.casefold() in {a.casefold() for a in aliases}]
    if len(candidates)!=1: raise ValueError(f'Expected one site directory for {aliases}, found {len(candidates)}')
    groups={}
    for p in candidates[0].rglob('*'):
        if p.is_file() and p.suffix.lower() in ('.tif','.tiff') and p.stem.isdigit():
            groups.setdefault(p.parent,[]).append(p)
    if len(groups)!=1: raise ValueError(f'Expected one TIFF folder, found {len(groups)}')
    folder,files=next(iter(groups.items()))
    numbers=[int(p.stem) for p in files]
    if len(numbers)!=len(set(numbers)) or set(numbers) not in (set(range(1,13)),set(range(1,5))):
        raise ValueError('Expected unique numbered set 1..12 or 1..4')
    # Downstream scripts currently use exact n.tif naming: do not rename raw data.
    if any(p.name!=f'{int(p.stem)}.tif' for p in files):
        raise ValueError('Requires exact n.tif names; raw data left unchanged')
    return folder,sorted(files,key=lambda p:int(p.stem))


def transferred_geometry(files,template):
    reference=next(r for r in template['records'] if Path(r['file']).name=='1.tif')
    f=next(f for f in reference['fits'] if f['threshold_fraction']==.05)
    h0,w0=reference['shape'][:2]; flags=[];records=[]
    for path in files:
        with tifffile.TiffFile(path) as t: shape=list(t.pages[0].shape)
        h,w=shape[:2]
        if len(shape)!=3 or shape[-1]!=3: raise ValueError('Expected RGB TIFF')
        scale=min(h,w)/min(h0,w0)
        # Centre fractions describe image placement; radius scales with short side.
        # Unknown crop changes are NOT corrected silently.
        fit=dict(threshold_fraction=.05,center_x_px=f['center_x_px']/w0*w,
            center_y_px=f['center_y_px']/h0*h,rim_radius_px=f['rim_radius_px']*scale)
        fit['r60_equisolid_nominal_px']=fit['rim_radius_px']/np.sqrt(2)
        if abs(w/h-w0/h0)>.02:
            flags.append('different aspect ratio: transferred geometry requires explicit crop review')
        records.append(dict(file=str(path.resolve()),shape=shape,sha256=digest(path),fits=[fit]))
    return dict(status='nominal_geometry_transfer_requires_review',channel=template['channel'],
        source='Cottage image 1 accepted 5% fit; short-side scale; fractional centre',
        records=records,raw_images_modified=False),sorted(set(flags))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--geometry-reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();source=a.root.resolve();out=a.output.resolve()
    if source==out or source in out.parents or out in source.parents: p.error('Output must be separate from raw root')
    if out.exists() and not a.resume: p.error('Output exists: use --resume to reuse verified completed sites')
    template=json.loads(a.geometry_reference.read_text())
    refs={r['site_id']:r for r in csv.DictReader((ROOT/'data/reference/prc_spider_18_sites_verified.csv').open())}
    out.mkdir(parents=True,exist_ok=True)
    rows=[];status=[]
    for name,(site_id,aliases) in SITES.items():
        print(f'[{len(status)+1}/18] {name}',flush=True)
        site=out/name;site.mkdir(exist_ok=True)
        try:
            folder,files=locate(source,aliases)
            file_hashes={str(f.resolve()):digest(f) for f in files}
            fingerprint=hashlib.sha256(json.dumps(dict(inputs=file_hashes,template=digest(a.geometry_reference),
                scripts={n:digest(ROOT/'src'/n) for n in ['batch_spider_prc_sites.py','analyse_spider_sky_rotation.py',
                    'audit_spider_polarization.py','reconstruct_spider_eye_contrast.py']}),sort_keys=True).encode()).hexdigest()
            marker=site/'complete.json'
            if a.resume and marker.exists() and json.loads(marker.read_text())['fingerprint']==fingerprint:
                saved=json.loads(marker.read_text());rows.extend(saved['rows']);status.append(saved['status'])
                print('  Reused completed site (inputs and scripts verified)',flush=True)
            else:
                if (site/'registration').exists() or (site/'polarization').exists():
                    raise ValueError('Partial or changed run exists: use a fresh batch output directory for this site')
                if name=='Cottage':
                    geometry=template;flags=[]
                    for r in geometry['records']:
                        if file_hashes.get(str((folder/Path(r['file']).name).resolve()))!=r['sha256']:
                            raise ValueError('Cottage template/input mismatch')
                else:
                    geometry,flags=transferred_geometry(files,template)
                    flags.append('geometry transferred from Cottage: inspect rim and field overlay')
                geometry_path=site/'geometry.json';geometry_path.write_text(json.dumps(geometry,indent=2)+'\n')
                commands=[['analyse_spider_sky_rotation.py','--geometry',str(geometry_path),'--folder',str(folder),
                    '--site-id',site_id,'--output',str(site/'registration')],
                    ['audit_spider_polarization.py','--run',str(site/'registration'),
                    '--angular-radius-deg','30','--output',str(site/'polarization')]]
                for command in commands:
                    with (site/(command[0]+'.log')).open('w') as log:
                        subprocess.run([sys.executable,str(ROOT/'src'/command[0]),*command[1:]],stdout=log,stderr=subprocess.STDOUT,check=True)
                rotation=json.loads((site/'registration/rotation_analysis.json').read_text())
                flagged=sum(bool(r['flags']) for r in rotation['rotations'])
                if flagged: flags.append(f'{flagged} exposures have registration flags')
                summary=json.loads((site/'polarization/polarization_summary.json').read_text())
                site_rows=[]
                for method,values in summary['methods'].items():
                    for mode,result in values['fixed_eye_results'].items():
                        reference=float(refs[site_id]['prc_sigg']);value=result['prc']
                        row=dict(site=name,site_id=site_id,replicates=len(files)//4,method=method,response_mode=mode,
                            full_field_deg=60,angular_radius_deg=30,radius_px=summary['radius_px'],
                            prc=value,reference_prc=reference,difference=value-reference,
                            percent_difference=100*(value-reference)/reference,
                            mean_dop=values['physical_mean_dop'],excluded_fraction=result['excluded_fraction'],
                            flagged_exposures=flagged,status='RECONSTRUCTED_REQUIRES_REVIEW',flags='; '.join(flags))
                        site_rows.append(row)
                site_status=dict(site=name,status='processed_requires_review',flags=flags)
                marker.write_text(json.dumps(dict(fingerprint=fingerprint,rows=site_rows,status=site_status),indent=2)+'\n')
                rows.extend(site_rows);status.append(site_status)
                print('  Processed; diagnostics retained',flush=True)
        except (ValueError,OSError,KeyError,StopIteration,subprocess.CalledProcessError) as e:
            status.append(dict(site=name,status='failed_or_missing',error=str(e)))
            print(f'  FAILED: {e}',flush=True)
        # Incremental checkpoints; failure at one site never hides other sites.
        (out/'batch_status.json').write_text(json.dumps(dict(sites=status,raw_images_modified=False,
            historical_fits_changed=False,full_field_deg=60,angular_radius_deg=30,
            excluded_extra_folders=['bristol','sandhammeren','wells']),indent=2)+'\n')
        if rows:
            with (out/'all_sites_prc.csv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(f'Finished: inspect {out}/batch_status.json and all_sites_prc.csv',flush=True)
    if any(s['status']=='failed_or_missing' for s in status): sys.exit(1)

if __name__=='__main__':main()
