"""Re-warp locked linear-green exposures; audit Stokes and export camera maps.

No landmark refit, per-frame normalization, celestial rotation or sky-model fit.
Existing locks and products are read only. Output requires a fresh directory.
"""
import argparse
import copy
import json
import shutil
from pathlib import Path
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reconstruct_spider_eye_contrast import digest, tiff_channel
from lock_spider_site_registration import warp
from correct_spider_analyzer_sign import evidence


def quartet_stokes(frames, sign):
    z, f, n, o = frames
    return z+n, .5*(z+f+n+o), z-n, sign*(f-o)


def stats(v, mask):
    v = v[mask & np.isfinite(v)]
    return dict(count=int(v.size), median=float(np.median(v)) if v.size else None,
                p95_absolute=float(np.percentile(abs(v), 95)) if v.size else None)


def run(lock_path, site, out, raw_manifest_root=None):
    base = lock_path.parent
    checks = json.loads((base/'reference_sha256.json').read_text())
    if checks.get(lock_path.name) != digest(lock_path):
        raise ValueError('Lock checksum mismatch')
    lock = json.loads(lock_path.read_text())['sites'][site]
    if lock.get('radiometry_status') != 'black_subtracted_linear_camera_green':
        raise ValueError('Requires registered linear camera-green RAW extraction')
    raw_audit = None
    if raw_manifest_root is not None:
        matches = [p for p in raw_manifest_root.rglob('linear_raw_manifest.json')
                   if digest(p) == lock.get('linear_raw_manifest_sha256')]
        if not matches:
            raise ValueError('Matching RAW preparation manifest not found')
        manifest = json.loads(matches[0].read_text())
        from prepare_spider_linear_raw import read_metadata
        for frame in manifest['frames']:
            raw = Path(frame['source_path'])
            if digest(raw) != frame['source_sha256']:
                raise ValueError('Original RAW bytes changed: '+str(raw))
            meta = read_metadata(raw)
            for key in ('datetime_original', 'exposure_seconds', 'iso', 'f_number'):
                if meta[key] != frame[key]:
                    raise ValueError('Original RAW metadata changed: '+str(raw))
        raw_audit = dict(manifest=str(matches[0]), sha256=digest(matches[0]),
                         original_raw_hashes_and_metadata_verified=True,
                         preparation=manifest['processing'], frame_count=len(manifest['frames']))
    saved = lock.get('analyzer_stokes_convention')
    unsigned = copy.deepcopy(lock)
    unsigned.pop('analyzer_stokes_convention', None)
    sign, trace = evidence(unsigned)
    if saved and saved['nominal_sequence_sign'] != sign:
        raise ValueError('Saved analyzer sign disagrees with landmark evidence')
    transforms = {int(Path(t['filename']).stem):t for t in lock['transforms']}
    records = {int(Path(t['path']).stem):t for t in lock['raw_inputs']}
    groups = lock['analyzer_groups']
    numbers = sorted({n for g in groups for n in g})
    c = np.asarray(lock['reference_center_xy_px'])
    r = lock['rim_radius_px']
    sums = None
    diagnostics = []
    hashes = {}
    for group in groups:
        frames = []
        for n in group:
            rec = records[n]
            path = Path(rec['path'])
            if digest(path) != rec['sha256']:
                raise ValueError('Prepared RAW image checksum changed: '+str(path))
            hashes[str(path)] = rec['sha256']
            im, _ = tiff_channel(path, None)
            tr = transforms[n]
            mapping = tr.get('fisheye_mapping') or np.asarray(tr['source_to_reference_xy'])
            frames.append(im if n == 1 else warp(im, mapping))
        ipair, iall, q, u = quartet_stokes(frames, sign)
        y, x = np.indices(ipair.shape)
        field = np.hypot(x-c[0], y-c[1]) < lock['eye_radius_px']
        with np.errstate(divide='ignore', invalid='ignore'):
            closure = ((frames[0]+frames[2])-(frames[1]+frames[3]))/(2*iall)
        diagnostics.append(dict(exposures=group, pair_sum_closure=stats(closure, field)))
        if sums is None:
            sums = [np.zeros_like(ipair, dtype=float) for _ in range(4)]
        for total, v in zip(sums, (ipair, iall, q, u)):
            total += v
    ipair, iall, q, u = [v/len(groups) for v in sums]
    oldpath = base/site/'registration/linear_stokes.npz'
    if checks.get(str(oldpath.relative_to(base))) != digest(oldpath):
        raise ValueError('Saved Stokes checksum mismatch')
    with np.load(oldpath, allow_pickle=False) as z:
        for k, value in [('I', ipair), ('Q', q), ('U', u if saved else sign*u)]:
            if not np.allclose(z[k], value, rtol=1e-10, atol=1e-9, equal_nan=True):
                raise ValueError('Recomputation disagrees with saved '+k)
    out.mkdir(parents=True, exist_ok=False)
    disc = np.hypot(x-c[0], y-c[1]) < .985*r
    fig, axes = plt.subplots(2, 2, figsize=(9, 8), layout='constrained')
    summaries = {}
    arrays = dict(I_pair=ipair, I_four_angle=iall, Q=q, U=u)
    for row, (name, intensity) in enumerate([('pair', ipair), ('four_angle', iall)]):
        with np.errstate(divide='ignore', invalid='ignore'):
            d = np.hypot(q, u)/intensity
        angle = np.rad2deg(.5*np.arctan2(u, q))
        valid = disc & np.isfinite(d) & (intensity > 0) & (d <= 1)
        avalid = valid & (np.hypot(q, u) > 0)
        arrays['DoLP_'+name] = d
        arrays['AoP_deg_'+name] = angle
        arrays['valid_'+name] = valid
        for v, label, mask, col in [(d, 'DoLP', valid, 0), (angle, 'AoP_deg', avalid, 1)]:
            v = np.where(mask, v, np.nan)
            tifffile.imwrite(out/(label+'_'+name+'.tif'), v.astype(np.float32), metadata=None)
            ax = axes[row, col]
            im = ax.imshow(v, cmap='viridis' if col == 0 else 'twilight', vmin=0 if col == 0 else -90, vmax=1 if col == 0 else 90)
            ax.set_title(name+': '+label)
            ax.axis('off')
            fig.colorbar(im, ax=ax, shrink=.7)
        summaries[name] = dict(central_dolp=stats(d, field), nonphysical_central_pixels=int((field & ~valid).sum()))
    np.savez_compressed(out/'measured_stokes_and_polarization.npz', **arrays)
    fig.suptitle(site+'; camera frame; analyzer zero unknown; no model fit')
    fig.savefig(out/'DoLP_AoP_camera.png', dpi=180)
    plt.close(fig)
    edge = base/site/'registration_edges.png'
    if edge.exists():
        shutil.copyfile(edge, out/'saved_registration_edges.png')
    report = dict(site=site, source_lock_sha256=digest(lock_path), source_stokes_sha256=digest(oldpath), script_sha256=digest(__file__), input_sha256=hashes,
        original_raw_audit=raw_audit, nominal_sequence_sign=sign, landmark_angle_evidence=trace, quartets=diagnostics, recomputation_matches_saved_stokes=True,
        summaries=summaries, equations='analyzer measurement J(theta)=(I+Q cos(2theta)+U sin(2theta))/2; pair I=J0+J90; four-angle I=sum(J)/2; Q=J0-J90; U=sign*(J45-J135); DoLP=hypot(Q,U)/I; AoP=atan2(U,Q)/2',
        limitations=['Signed nominal analyzer angles; exact roll departures remain uncorrected.', 'No fitted frame gains or offsets. Pair closure tests ideal-model consistency, not its cause.', 'Four-angle equal-weight least squares is a diagnostic alternative; not an automatic correction of inconsistent acquisition.', 'Replicate I,Q,U averaged before DoLP/AoP; unsupported pixels remain invalid.', 'Camera plane only, unknown analyzer zero; no north transform or calibrated off-axis Mueller transport.', 'Saved alignment overlays copied for reference; no new visual approval inferred.', 'Original RAW bytes and metadata verified only when original_raw_audit is present.'])
    (out/'stokes_audit.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(site+': recomputation matches; exported '+str(out), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock', type=Path, required=True)
    p.add_argument('--sites', nargs='+', default=['Cottage', 'museum'])
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--raw-manifest-root', type=Path, help='Find matching preparation manifests and verify original NEF bytes/metadata')
    a = p.parse_args()
    if a.output.exists():
        p.error('Use a fresh output directory')
    a.output.mkdir(parents=True)
    result = {}
    for site in a.sites:
        try:
            run(a.reference_lock, site, a.output/site, a.raw_manifest_root)
            result[site] = dict(status='exported_not_new_visual_approval')
        except (OSError, ValueError, KeyError) as exc:
            result[site] = dict(status='pending', reason=str(exc))
            print(site+': PENDING '+str(exc), flush=True)
    (a.output/'audit_status.json').write_text(json.dumps(result, indent=2)+'\n')
    if any(r['status']=='pending' for r in result.values()):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
