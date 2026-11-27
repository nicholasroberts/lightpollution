"""Inspect ChurchHill/zoology saved registrations before natural-pattern fitting."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from correct_spider_analyzer_sign import evidence

SITES=('ChurchHill','zoology')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def inspect(path,site):
    lock=json.loads(path.read_text())['sites'][site]
    checks=json.loads((path.parent/'reference_sha256.json').read_text())
    stokes=path.parent/site/'registration/linear_stokes.npz'
    for source in (path,stokes):
        if checks.get(str(source.relative_to(path.parent)))!=digest(source):
            raise ValueError('Checksum differs: '+str(source))
    convention=lock.get('analyzer_stokes_convention')
    if convention:
        sign=convention['nominal_sequence_sign']
        trace=convention.get('quartet_rotation_evidence',[])
        sign_status='already_applied_do_not_repeat'
    else:
        sign,trace=evidence(lock)
        sign_status='inferred_from_landmarks_not_yet_applied'
    anchor=lock.get('astronomical_anchor')
    row=dict(site=site,lock_path=str(path),lock_sha256=digest(path),
        stokes_path=str(stokes),stokes_sha256=digest(stokes),
        status=lock.get('status'),radiometry_status=lock.get('radiometry_status'),
        moon_pixel_xy=lock.get('moon_pixel_xy'),astronomical_anchor=anchor,
        moon_position_status=lock.get('moon_position_status'),
        moon_localization=lock.get('moon_localization'),
        visual_alignment_status=lock.get('visual_alignment_status'),
        north_orientation_status=lock.get('north_orientation_status'),
        celestial_orientation=lock.get('celestial_orientation'),
        analyzer_sequence_sign=sign,analyzer_sign_status=sign_status,
        analyzer_rotation_evidence=trace,
        polarization_fitted=False,registration_changed=False)
    with np.load(stokes,allow_pickle=False) as z:
        row['stokes_shape']=list(z['I'].shape)
    raw=next((r for r in lock.get('raw_inputs',[]) if Path(r['path']).stem=='1'),None)
    row['reference_image']=raw
    return row

def run(a):
    if a.output.exists():raise ValueError('Use a fresh output directory')
    rows=[]
    # Report all candidates; never select an arbitrary newer lock.
    for site in SITES:
        candidates=[]
        for path in sorted(a.processed.rglob('locked_site_references.json')):
            try:
                payload=json.loads(path.read_text())
                if site not in payload.get('sites',{}):continue
                result=inspect(path,site)
                candidates.append(result)
            except (OSError,ValueError,KeyError) as exc:
                candidates.append(dict(site=site,lock_path=str(path),inspection_error=str(exc)))
        rows.append(dict(site=site,candidates=candidates))
    a.output.mkdir(parents=True)
    (a.output/'registration_readiness.json').write_text(json.dumps(rows,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import tifffile
    drawable=[]
    seen=set()
    for row in rows:
        for candidate in row['candidates']:
            ref=candidate.get('reference_image')
            if ref and (row['site'],ref['sha256']) not in seen:
                seen.add((row['site'],ref['sha256']));drawable.append(candidate)
    if drawable:
        fig,axes=plt.subplots(len(drawable),1,figsize=(12,7*len(drawable)),squeeze=False,layout='constrained')
        for ax,row in zip(axes[:,0],drawable):
            ref=row['reference_image'];path=Path(ref['path'])
            try:
                if digest(path)!=ref['sha256']:raise ValueError('Reference image checksum differs')
                image=tifffile.imread(path)
                lo,hi=np.nanpercentile(image,[5,99.9])
                image=np.log1p(100*np.clip((image-lo)/max(hi-lo,1e-12),0,1))
                ax.imshow(image,cmap='gray')
                moon=row.get('moon_pixel_xy')
                if moon is not None:ax.plot(*moon,'r+',ms=14)
                ax.set_title(row['site']+' — saved candidate Moon; identity requires review')
            except (OSError,ValueError) as exc:ax.text(.1,.5,str(exc),transform=ax.transAxes)
            ax.set_xlabel('Pixel x');ax.set_ylabel('Pixel y')
        fig.savefig(a.output/'reference_moon_check.png',dpi=180);plt.close(fig)
    print(json.dumps(rows,indent=2))
    print('Saved',a.output)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--processed',type=Path,default=Path('data/processed'))
    p.add_argument('--output',type=Path,default=Path('data/processed/paper1_supp03_readiness_01'))
    a=p.parse_args()
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))

if __name__=='__main__':main()
