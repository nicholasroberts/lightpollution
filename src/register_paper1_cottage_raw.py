"""Fresh Cottage green-RAW registration for camera-frame 1A review; ephemeris pending."""
import argparse,json,subprocess,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--linear-raw-manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--registration-model',choices=['fisheye-shared','fisheye'],default='fisheye-shared',
        help='Existing constrained geometric model; validation thresholds are unchanged')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    manifest=json.loads(a.linear_raw_manifest.read_text())
    if manifest.get('camera_channel')!='green' or len(manifest.get('frames',[]))!=12:
        p.error('Requires all twelve prepared green frames')
    if manifest.get('analyzer_groups')!=[[1,4,7,10],[2,5,8,11],[3,6,9,12]]:
        p.error('Expected the three investigator-confirmed quartets')
    # Registration uses the previous report only for starting reference geometry.
    # Do not inherit the old TIFF capture time, Moon pixel or celestial rotation.
    a.output.mkdir(parents=True)
    pending=a.output/'no_inherited_ephemeris'/'hps'
    pending.mkdir(parents=True)
    (pending/'provenance.json').write_text(json.dumps(dict(
        sites=[dict(site='Cottage',anchor=None)],
        note='Camera-frame registration only; new NEF capture time and Moon pixel require review. No inherited TIFF anchor.'),
        indent=2)+'\n')
    subprocess.run([sys.executable,str(Path(__file__).with_name('lock_spider_site_registration.py')),
        '--batch',str(a.batch),'--model-run',str(pending.parent),
        '--sites','Cottage','--linear-raw-manifest',str(a.linear_raw_manifest),
        '--raw-time-mapping',str(Path(__file__).resolve().parents[1]/'config/paper1/cottage_raw_mapping.json'),
        '--registration-model',a.registration_model,'--output',str(a.output/'registration')],check=True)
    registered=json.loads((a.output/'registration/locked_site_references.json').read_text())['sites']['Cottage']
    if registered.get('status')!='numerically_locked_awaiting_visual_review':
        p.error('Cottage registration remains pending: '+registered.get('reason',registered.get('status','unknown'))
                +'. Inspect registration outputs; Stokes audit has not been attempted.')
    subprocess.run([sys.executable,str(Path(__file__).with_name('audit_spider_measured_stokes.py')),
        '--reference-lock',str(a.output/'registration/locked_site_references.json'),
        '--sites','Cottage','--raw-manifest-root',str(a.linear_raw_manifest.parent),
        '--output',str(a.output/'audit')],check=True)
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cmcrameri import cm
    lock=json.loads((a.output/'registration/locked_site_references.json').read_text())['sites']['Cottage']
    with np.load(a.output/'audit/Cottage/measured_stokes_and_polarization.npz',allow_pickle=False) as z:
        intensity,q,u=[np.array(z[k]) for k in ['I_pair','Q','U']]
    with np.errstate(divide='ignore',invalid='ignore'):d=np.hypot(q,u)/intensity
    angle=np.rad2deg(.5*np.arctan2(u,q))
    y,x=np.indices(d.shape);c=np.asarray(lock['reference_center_xy_px'])
    field=np.hypot(x-c[0],y-c[1])<lock['eye_radius_px']
    disc=np.hypot(x-c[0],y-c[1])<.985*lock['rim_radius_px']
    valid=(intensity>0)&np.isfinite(d)&(d>=0)&(d<=1)
    values=d[field&valid]
    summary=dict(site='Cottage',mean_dolp=float(values.mean()),spatial_sd=float(values.std(ddof=1)),
        included_pixels=int(values.size),excluded_field_pixels=int((field&~valid).sum()),
        definition='Mean per-pixel camera-green DoLP and sample spatial SD in zenith-centred 60-degree full field; pair estimator; before north resampling',
        orientation='camera frame; north and Moon anchor pending',all_twelve_frames_used=True,
        source_lock_sha256=__import__('hashlib').sha256((a.output/'registration/locked_site_references.json').read_bytes()).hexdigest())
    fig,axes=plt.subplots(1,2,figsize=(10,5.3),layout='constrained')
    for ax,v,mask,cmap,lo,hi,label in [
        (axes[0],d,valid,plt.get_cmap('viridis'),0,1,'Measured green DoLP'),
        (axes[1],angle,valid&(np.hypot(q,u)>0),cm.romaO,-90,90,'Measured green AoP')]:
        im=ax.imshow(np.where(disc&mask,v,np.nan),cmap=cmap,vmin=lo,vmax=hi)
        ax.add_patch(plt.Circle(c,lock['eye_radius_px'],fill=False,color='white'))
        ax.set_title(label);ax.axis('off');fig.colorbar(im,ax=ax,shrink=.7)
    fig.suptitle('Cottage — camera frame; north/Moon check pending\n'
        +format(summary['mean_dolp'],'.3f')+' ± '+format(summary['spatial_sd'],'.3f')+' (1 spatial SD, full 60° field)')
    fig.savefig(a.output/'Cottage_camera_DoLP_AoP.png',dpi=200);plt.close(fig)
    (a.output/'Cottage_camera_field_statistics.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Review '+str(a.output/'registration/Cottage/registration_edges.png'))
    print('Camera DoLP review: '+str(a.output/'registration/Cottage/registered_dop_and_north.png'))
    print('North remains unset until a fresh Moon/ephemeris check.')


if __name__=='__main__':main()
