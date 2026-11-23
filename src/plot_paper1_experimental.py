"""Figure 1A: experimental green DoLP/AoP from locked Stokes and saved geometry.

No registration, parity, analyzer-zero or scientific parameter is fitted here.
Each site needs its own saved fit report to specify the reviewed display parity.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def block(v, n):
    h, w = v.shape
    return v[:h//n*n, :w//n*n].reshape(h//n, n, w//n, n).mean((1, 3))


def transform_stokes(q, u, a):
    """Transform polarization tensor with the same orthogonal image transform."""
    t = np.empty(q.shape + (2, 2))
    t[..., 0, 0] = q
    t[..., 1, 1] = -q
    t[..., 0, 1] = t[..., 1, 0] = u
    out = np.einsum('ij,...jk,lk->...il', a, t, a)
    return out[..., 0, 0], out[..., 0, 1]


def prepare(lock_path, fit_path, site, n):
    lock = json.loads(lock_path.read_text())['sites'][site]
    fit = json.loads(fit_path.read_text())
    sp = lock_path.parent/site/'registration/linear_stokes.npz'
    checks = json.loads((lock_path.parent/'reference_sha256.json').read_text())
    for p in (lock_path, sp):
        if checks.get(str(p.relative_to(lock_path.parent))) != sha(p):
            raise ValueError('Lock checksum mismatch: '+str(p))
    if fit['site'] != site or fit['source_lock_sha256'] != sha(lock_path) or fit['stokes_sha256'] != sha(sp):
        raise ValueError('Fit geometry and locked Stokes differ: '+site)
    if lock.get('camera_channel') != 'green' or lock.get('radiometry_status') != 'black_subtracted_linear_camera_green':
        raise ValueError('Current linear green reconstruction required: '+site)
    if lock.get('moon_pixel_xy') is None:
        raise ValueError('Observed Moon anchor required: '+site)
    with np.load(sp, allow_pickle=False) as z:
        i, q, u = [block(z[k], n) for k in ('I', 'Q', 'U')]
    c = np.asarray(lock['reference_center_xy_px'])
    r = lock['rim_radius_px']
    m = np.asarray(fit['display_camera_to_north_xy'])
    if m.shape != (2, 3) or not np.allclose(m[:, :2]@m[:, :2].T, np.eye(2), atol=1e-6):
        raise ValueError('Expected orthogonal saved display transform')
    # Saved display is north-up/E-right in down-positive pixel coordinates.
    # Reverse both display axes to obtain north-up/E-left Cartesian sky view.
    a = -m[:, :2]
    y, x = np.indices(i.shape)
    xy = np.stack([x*n+(n-1)/2, y*n+(n-1)/2], -1)
    pos = -(xy@m[:, :2].T+m[:, 2]-c)/r
    qd, ud = transform_stokes(q, u, a)
    with np.errstate(divide='ignore', invalid='ignore'):
        d = np.hypot(q, u)/i
    valid = (np.linalg.norm(xy-c, axis=-1) < .985*r)&(i > 0)&np.isfinite(d)&(d <= 1)
    angle = np.rad2deg(.5*np.arctan2(ud, qd))
    moon = -(m[:, :2]@np.asarray(lock['moon_pixel_xy'])+m[:, 2]-c)/r
    return pos, np.where(valid, d, np.nan), np.where(valid&(d >= .05), angle, np.nan), moon, lock['eye_radius_px']/r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock', type=Path, required=True)
    p.add_argument('--fit-root', type=Path, required=True, help='Contains SITE/fit/twilight_moon_fit.json')
    p.add_argument('--sites', nargs='+', default=['Cottage', 'museum'])
    p.add_argument('--bin-size', type=int, default=8)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.bin_size < 1 or a.output.exists():
        p.error('Positive bin size and fresh output directory required')
    try:
        paths = [a.fit_root/s/'fit/twilight_moon_fit.json' for s in a.sites]
        panels = [prepare(a.reference_lock, f, s, a.bin_size) for s, f in zip(a.sites, paths)]
    except (OSError, KeyError, ValueError) as exc:
        p.error(str(exc))
    plt.rcParams.update({'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig, axes = plt.subplots(len(panels), 2, squeeze=False, figsize=(7, 3.6*len(panels)), layout='constrained')
    for row, (site, panel) in enumerate(zip(a.sites, panels)):
        pos, d, angle, moon, eye = panel
        for col, (v, cmap, lo, hi, title) in enumerate([(d, 'viridis', 0, 1, 'Measured DoLP'), (angle, 'twilight', -90, 90, 'Measured AoP; unknown constant zero')]):
            ax = axes[row, col]
            im = ax.pcolormesh(pos[..., 0], pos[..., 1], v, shading='nearest', cmap=cmap, vmin=lo, vmax=hi, rasterized=True)
            ax.add_patch(Circle((0, 0), eye, fill=False, color='white', lw=.8))
            ax.plot(*moon, 'kx', ms=6)
            for x, y, text in [(0, 1.06, 'N'), (-1.06, 0, 'E'), (0, -1.06, 'S'), (1.06, 0, 'W')]:
                ax.text(x, y, text, ha='center', va='center', fontsize=9)
            ax.set(aspect='equal', xlim=(-1.15, 1.15), ylim=(-1.15, 1.15), title=site+'\n'+title)
            ax.axis('off')
            fig.colorbar(im, ax=ax, shrink=.75, label='DoLP' if col == 0 else 'Degrees modulo 180')
    a.output.mkdir(parents=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(a.output/('fig01a_experimental.'+ext), dpi=300)
    plt.close(fig)
    provenance = dict(sites=a.sites, source_lock_sha256=sha(a.reference_lock), fit_reports_sha256={s:sha(f) for s,f in zip(a.sites, paths)}, script_sha256=sha(__file__), bin_size=a.bin_size,
        interpretation='Experimental camera-green only. Full disc includes foreground; no automatic clear-sky classification. Circle is central 60 degree full field. AoP shown only for physical DoLP >=0.05; constant analyzer zero unknown. Camera-plane tensor transformed with image; off-axis transport uncalibrated. No fit performed.')
    (a.output/'fig01a_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
