"""Correct nominal analyzer U sign from independently saved camera rotations.

Assumes the rear analyzer is fixed to and rotates with the camera/lens.
Uses source-to-reference landmark rotations, never agreement with a sky model.
"""
import argparse
import copy
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def roll(transform):
    if 'fisheye_mapping' in transform:
        r=np.asarray(transform['fisheye_mapping']['ray_rotation_source_to_reference'])
    else:r=np.asarray(transform['source_to_reference_xy'])[:,:2]
    return float(np.rad2deg(np.arctan2(r[1,0],r[0,0])))


def evidence(lock):
    if lock.get('analyzer_stokes_convention'):raise ValueError('Already has an explicit Stokes convention; refusing a second correction')
    if lock.get('status')!='numerically_locked_awaiting_visual_review':raise ValueError('Numerical registration pending')
    transforms={int(Path(t['filename']).stem):t for t in lock['transforms']}
    records=[];signs=[]
    nominal=np.array([0.,45.,90.,135.])
    for group in lock['analyzer_groups']:
        if len(group)!=4:raise ValueError('Require explicit four-angle quartets')
        if not all(transforms[n]['quality']['numerical_pass'] for n in group):raise ValueError('Nonpassing exposure in quartet')
        actual=np.array([roll(transforms[n]) for n in group]);actual-=actual[0]
        errors={sign:(actual-sign*nominal+90)%180-90 for sign in [1,-1]}
        acceptable=[sign for sign,v in errors.items() if np.max(abs(v))<10]
        if len(acceptable)!=1:raise ValueError('Camera roll does not unambiguously establish signed nominal sequence')
        sign=acceptable[0];signs.append(sign)
        records.append(dict(exposures=group,source_to_reference_roll_relative_deg=actual.tolist(),
            signed_nominal_angles_deg=(sign*nominal).tolist(),
            maximum_nominal_angle_deviation_deg=float(np.max(abs(errors[sign]))),
            positive_sequence_maximum_deviation_deg=float(np.max(abs(errors[1]))),
            negative_sequence_maximum_deviation_deg=float(np.max(abs(errors[-1])))))
    if not signs or len(set(signs))!=1:raise ValueError('Quartets disagree about signed analyzer sequence')
    return signs[0],records


def apply(a):
    source=a.reference_lock.parent
    checks=json.loads((source/'reference_sha256.json').read_text())
    for rel,sha in checks.items():
        path=source/rel
        if not path.resolve().is_relative_to(source.resolve()) or digest(path)!=sha:
            raise ValueError('Source checksum mismatch: '+str(path))
    if checks.get(a.reference_lock.name)!=digest(a.reference_lock):raise ValueError('Source lock not in checksum manifest')
    original=json.loads(a.reference_lock.read_text());prepared=[]
    for site in a.sites:
        lock=copy.deepcopy(original['sites'][site]);sign,trace=evidence(lock)
        path=source/site/'registration/linear_stokes.npz'
        with np.load(path) as z:arrays={k:np.array(z[k]) for k in z.files}
        if any(k not in arrays for k in ['I','Q','U']):raise ValueError('Missing Stokes component')
        if any(arrays[k].shape!=arrays['I'].shape for k in ['Q','U']):raise ValueError('Stokes shapes differ')
        prepared.append((site,lock,sign,trace,arrays,path))
    a.output.mkdir(parents=True,exist_ok=False);locks={};summary=[]
    for site,lock,sign,trace,arrays,path in prepared:
        shutil.copytree(source/site,a.output/site)
        before={k:v.copy() for k,v in arrays.items()};arrays['U']=sign*arrays['U']
        out=a.output/site/'registration/linear_stokes.npz';np.savez_compressed(out,**arrays)
        with np.load(out) as z:
            for key in ['I','Q']:np.testing.assert_array_equal(z[key],before[key])
            np.testing.assert_array_equal(z['U'],sign*before['U'])
            np.testing.assert_array_equal(np.hypot(z['Q'],z['U']),np.hypot(before['Q'],before['U']))
        convention=dict(basis='reference-image analyzer plane: x right, y down',
            analyzer_mount_assumption='rear analyzer fixed to and co-rotating with camera/lens',
            decision_source='signed source-to-reference landmark rotations; no polarization or astronomical objective',
            nominal_sequence_sign=sign,applied_U_multiplier=sign,source_stokes_sha256=digest(path),
            source_lock_sha256=digest(a.reference_lock),quartet_rotation_evidence=trace,
            I_Q_and_DoLP_unchanged=True,exposure_alignment_unchanged=True,exact_roll_angle_refit=False,
            analyzer_zero_status='unknown constant relative to reference-image x axis',
            limitations=['Signed nominal angles only; measured departures from 45-degree steps are not corrected.',
                'Near-axis rotation evidence does not provide full off-axis lens/analyzer Mueller calibration.',
                'Camera-to-sky parity and publication-view reflection remain separate from analyzer sequence sign.'])
        lock['analyzer_stokes_convention']=convention
        if lock['transforms']!=original['sites'][site]['transforms']:raise ValueError('Exposure transforms changed')
        (a.output/site/'reference_coords.json').write_text(json.dumps(lock,indent=2,allow_nan=False)+'\n')
        reportpath=a.output/site/'registration/rotation_analysis.json'
        if reportpath.exists():
            report=json.loads(reportpath.read_text());report['analyzer_stokes_convention']=convention
            reportpath.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        locks[site]=lock;summary.append(dict(site=site,**convention))
    derived=dict(original,sites=locks,status='SIGNED_NOMINAL_ANALYZER_CONVENTION_FROM_LANDMARK_ROTATION',
        source_lock_sha256=digest(a.reference_lock),script_sha256=digest(__file__))
    (a.output/'locked_site_references.json').write_text(json.dumps(derived,indent=2,allow_nan=False)+'\n')
    (a.output/'analyzer_sign_summary.json').write_text(json.dumps(dict(sites=summary),indent=2)+'\n')
    checks={str(f.relative_to(a.output)):digest(f) for f in a.output.rglob('*') if f.is_file()}
    (a.output/'reference_sha256.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--sites',nargs='+',default=['WClawn','oxfordbussinesspark'])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try:apply(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))


if __name__=='__main__':main()
