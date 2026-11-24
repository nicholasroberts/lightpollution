"""Render measured subfigure 1A from saved north-up Stokes; no re-rotation."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cmcrameri import cm
from reconstruct_spider_eye_contrast import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--north',type=Path,required=True)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--site',default='museum')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    north_report=a.north/'north_rotation.json';audit_report=a.audit/'stokes_audit.json'
    nr=json.loads(north_report.read_text());ar=json.loads(audit_report.read_text())
    if nr['site']!=a.site or ar['site']!=a.site:
        p.error('Site mismatch')
    for report in [nr,ar]:
        if report['source_lock_sha256']!=digest(a.reference_lock):
            p.error('Source lock mismatch')
    if nr['source_audit_sha256']!=digest(audit_report):
        p.error('North display and audit report differ')
    camera_path=a.audit/'measured_stokes_and_polarization.npz'
    if nr['source_arrays_sha256']!=digest(camera_path):
        p.error('North display was generated from different measured arrays')
    lock=json.loads(a.reference_lock.read_text())['sites'][a.site]
    with np.load(a.north/'north_stokes.npz',allow_pickle=False) as z:
        d=np.array(z['DoLP']);angle=np.array(z['AoP_deg'])
        center=np.array(z['center_xy']);moon=np.array(z['moon_xy'])
    with np.load(camera_path,allow_pickle=False) as z:
        i,q,u=[z[k] for k in ['I_pair','Q','U']]
        with np.errstate(divide='ignore',invalid='ignore'):
            original=np.hypot(q,u)/i
    y,x=np.indices(i.shape);c=np.asarray(lock['reference_center_xy_px'])
    field=np.hypot(x-c[0],y-c[1])<lock['eye_radius_px']
    good=field&(i>0)&np.isfinite(original)&(original>=0)&(original<=1)
    if good.sum()<2:p.error('Insufficient valid camera-field pixels')
    mean=float(original[good].mean());sd=float(original[good].std(ddof=1))
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(10,5.3),layout='constrained')
    for col,(ax,v,cmap,lo,hi,label) in enumerate(zip(axes,[d,angle],[plt.get_cmap('viridis'),cm.romaO],[0,-90],[1,90],['Measured green DoLP','Measured green AoP'])):
        im=ax.imshow(v,cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
        ax.add_patch(plt.Circle(center,lock['eye_radius_px'],fill=False,color='white',lw=.9))
        ax.plot(*moon,'k+',ms=6)
        for tx,ty,t in [(.5,.99,'N'),(.01,.5,'E'),(.5,.01,'S'),(.99,.5,'W')]:
            ax.text(tx,ty,t,transform=ax.transAxes,ha='center',va='center')
        title=label
        if col==0:title+='\n'+f'{mean:.3f} ± {sd:.3f} (spatial SD)'
        ax.set_title(title,fontsize=11)
        r=lock['rim_radius_px'];ax.set_xlim(center[0]-1.07*r,center[0]+1.07*r);ax.set_ylim(center[1]+1.07*r,center[1]-1.07*r)
        ax.axis('off');fig.colorbar(im,ax=ax,shrink=.8,label='DoLP' if col==0 else 'Degrees modulo 180')
    fig.suptitle(a.site+' — measured polarization\nMoon-anchored north up, east left; absolute analyzer zero unknown',fontsize=12)
    fig.text(.015,.985,'A',ha='left',va='top',fontweight='bold',fontsize=18)
    a.output.mkdir(parents=True)
    for ext in ['png','pdf','svg']:
        fig.savefig(a.output/('fig01a_'+a.site+'.'+ext),dpi=300)
    plt.close(fig)
    report=dict(site=a.site,colormaps=dict(DoLP='viridis',AoP='romaO'),
        camera_field_mean_dolp=mean,camera_field_spatial_sd=sd,
        included_pixels=int(good.sum()),excluded_field_pixels=int((field&~good).sum()),
        statistic_definition='Unrotated camera-pixel DoLP mean and sample SD ddof=1; zenith-centred full 60-degree field; pair intensity; no added sky mask',
        source_lock_sha256=digest(a.reference_lock),source_audit_sha256=digest(audit_report),
        source_camera_arrays_sha256=digest(camera_path),source_north_arrays_sha256=digest(a.north/'north_stokes.npz'),
        source_north_report_sha256=digest(north_report),script_sha256=digest(__file__),
        numerical_processing_changed=False,new_rotation_or_resampling=False)
    (a.output/'figure_values.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'{a.site}: {mean:.6f} ± {sd:.6f}; saved {a.output}')


if __name__=='__main__':main()
