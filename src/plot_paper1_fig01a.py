"""Assemble Cottage and museum measured 1A from saved north-up arrays."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cmcrameri import cm
from paper1_figure_archive import archive_figure, digest


def load_site(site, north, audit, reference_lock):
    north, audit, reference_lock = map(Path, (north, audit, reference_lock))
    nr_path = north/'north_rotation.json'
    ar_path = audit/'stokes_audit.json'
    camera_path = audit/'measured_stokes_and_polarization.npz'
    nr, ar = [json.loads(p.read_text()) for p in (nr_path, ar_path)]
    for report in (nr, ar):
        if report['site'] != site or report['source_lock_sha256'] != digest(reference_lock):
            raise ValueError(site+': source site/lock mismatch')
    if nr['source_audit_sha256'] != digest(ar_path) or nr['source_arrays_sha256'] != digest(camera_path):
        raise ValueError(site+': north product differs from measured audit/arrays')
    lock = json.loads(reference_lock.read_text())['sites'][site]
    if lock.get('camera_channel') != 'green':
        raise ValueError(site+': camera-green lock required')
    with np.load(north/'north_stokes.npz', allow_pickle=False) as z:
        d, angle, center, moon = [np.array(z[k]) for k in ('DoLP','AoP_deg','center_xy','moon_xy')]
    with np.load(camera_path, allow_pickle=False) as z:
        i, q, u = [z[k] for k in ('I_pair','Q','U')]
        with np.errstate(divide='ignore', invalid='ignore'):
            original = np.hypot(q,u)/i
    y,x = np.indices(i.shape)
    c = np.asarray(lock['reference_center_xy_px'])
    field = np.hypot(x-c[0], y-c[1]) < lock['eye_radius_px']
    good = field & (i>0) & np.isfinite(original) & (original>=0) & (original<=1)
    if good.sum()<2:
        raise ValueError(site+': insufficient valid original camera-field pixels')
    values = dict(site=site, camera_field_mean_dolp=float(original[good].mean()),
        camera_field_spatial_sd=float(original[good].std(ddof=1)),
        included_pixels=int(good.sum()), excluded_field_pixels=int((field & ~good).sum()),
        astronomical_anchor=lock['astronomical_anchor'],
        sources={str(p):digest(p) for p in (reference_lock, ar_path, camera_path, nr_path, north/'north_stokes.npz')})
    return dict(d=d, angle=angle, center=center, moon=moon, lock=lock, values=values)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cottage-root', type=Path, default=Path('data/processed/paper1_cottage_north_01'))
    p.add_argument('--museum-north', type=Path, default=Path('data/processed/paper1_museum_north_01'))
    p.add_argument('--museum-audit', type=Path, default=Path('data/processed/paper1_museum_stokes_audit_01/museum'))
    p.add_argument('--museum-lock', type=Path, default=Path('data/processed/spider_registration_green_all_sites_01/runs/museum/locked_site_references.json'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--figure-copy', type=Path, required=True, help='Fresh version under figures/paper1; never overwritten')
    a = p.parse_args()
    if a.output.exists() or a.figure_copy.exists():
        p.error('Output and figure-copy must both be fresh versioned directories')
    try:
        sites = [load_site('Cottage', a.cottage_root/'north', a.cottage_root/'audit/Cottage',
                          a.cottage_root/'registration/locked_site_references.json'),
                 load_site('museum', a.museum_north, a.museum_audit, a.museum_lock)]
    except (OSError, ValueError, KeyError) as exc:
        p.error(str(exc))
    plt.rcParams.update({'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig, axes = plt.subplots(2,2,figsize=(10,10),layout='constrained')
    for row, item in enumerate(sites):
        vals, lock = item['values'], item['lock']
        for col, (key,cmap,lo,hi,label) in enumerate([
            ('d','viridis',0,1,'DoLP'), ('angle',cm.romaO,-90,90,'AoP')]):
            ax=axes[row,col]
            im=ax.imshow(item[key],cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
            ax.add_patch(plt.Circle(item['center'],lock['eye_radius_px'],fill=False,color='white',lw=.9))
            ax.plot(*item['moon'],'k+',ms=6)
            for tx,ty,t in [(.5,.99,'N'),(.01,.5,'E'),(.5,.01,'S'),(.99,.5,'W')]:
                ax.text(tx,ty,t,transform=ax.transAxes,ha='center',va='center')
            title=vals['site']+' — measured green '+label
            if col==0:
                title+='\n'+f"{vals['camera_field_mean_dolp']:.3f} ± {vals['camera_field_spatial_sd']:.3f} (spatial SD)"
            ax.set_title(title,fontsize=11)
            c=item['center'];r=lock['rim_radius_px']
            ax.set_xlim(c[0]-1.07*r,c[0]+1.07*r);ax.set_ylim(c[1]+1.07*r,c[1]-1.07*r)
            ax.axis('off')
            fig.colorbar(im,ax=ax,shrink=.8,label='DoLP' if col==0 else 'Degrees modulo 180')
    fig.suptitle('Measured polarization\nMoon-anchored north up, east left; absolute analyzer zeros unknown',fontsize=12)
    fig.text(.015,.985,'A',ha='left',va='top',fontweight='bold',fontsize=18)
    a.output.mkdir(parents=True)
    for ext in ('png','pdf','svg'):
        fig.savefig(a.output/('fig01a_Cottage_museum.'+ext),dpi=300)
    plt.close(fig)
    report=dict(subfigure='1A',rows=['Cottage','museum'],columns=['DoLP','AoP'],
        sites=[s['values'] for s in sites],colormaps=dict(DoLP='viridis',AoP='romaO'),
        dolp_colour_normalization='linear 0–1',aop_colour_limits_deg=[-90,90],
        statistic_definition='Unrotated camera-pixel DoLP mean and sample SD ddof=1; zenith-centred full 60-degree field; pair intensity; no added sky mask',
        script_sha256=digest(__file__),new_rotation_or_resampling=False,numerical_processing_changed=False)
    (a.output/'figure_values.json').write_text(json.dumps(report,indent=2)+'\n')
    with (a.output/'field_statistics.csv').open('w',newline='') as stream:
        keys=['site','camera_field_mean_dolp','camera_field_spatial_sd','included_pixels','excluded_field_pixels']
        writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader()
        writer.writerows({k:s['values'][k] for k in keys} for s in sites)
    archive_figure(a.output,a.figure_copy)
    print('Saved',a.output,'and verified figure snapshot',a.figure_copy)


if __name__=='__main__':
    main()
