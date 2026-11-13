"""Calculate Moon elevation/azimuth from investigator-confirmed local camera time.

Cottage Europe/Stockholm; Oxford Europe/London. Uses original-exposure EXIF
if available, otherwise explicitly reports TIFF DateTime fallback. Dates that
disagree with reference field date are flagged, not silently rewritten.
Predicted elevation constrains a ring in a level fisheye image; image north is
unknown, so a brightness candidate on that ring is still provisional.
"""
from pathlib import Path
from datetime import datetime,timezone
from zoneinfo import ZoneInfo
import argparse,csv,json
import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from wang2016_sky import lunar_track
from reconstruct_spider_eye_contrast import digest
ROOT=Path(__file__).resolve().parents[1]
SITES={'Cottage':('s14','Europe/Stockholm'),'oxfordbussinesspark':('s02','Europe/London'),'Skillinge':('s10','Europe/Stockholm'),'slabhouse':('s17','Europe/London')}


def local_to_utc(text,zone):
    parsed=datetime.strptime(text,'%Y:%m:%d %H:%M:%S').replace(tzinfo=ZoneInfo(zone))
    return parsed.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S'),parsed.utcoffset().total_seconds()/3600


def capture_time(path):
    import tifffile
    with tifffile.TiffFile(path) as tf:
        tags=tf.pages[0].tags
        exif=tags.get('ExifTag');values=exif.value if exif is not None else {}
        for key in ['DateTimeOriginal',36867,'DateTimeDigitized',36868]:
            if key in values:return str(values[key]),f'EXIF {key}'
        tag=tags.get('DateTime')
        if tag is None:raise ValueError('No camera timestamp found')
        return str(tag.value),'TIFF DateTime fallback: may be export metadata'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cottage-local-time',help='Explicit corrected local timestamp YYYY:MM:DD HH:MM:SS')
    p.add_argument('--oxford-local-time',help='Explicit corrected local timestamp YYYY:MM:DD HH:MM:SS')
    p.add_argument('--sites',nargs='+',choices=list(SITES),default=['Cottage','oxfordbussinesspark'])
    p.add_argument('--local-time-overrides',type=Path,help='JSON {SITE: YYYY:MM:DD HH:MM:SS}; use only independently confirmed dates/times')
    a=p.parse_args();out=a.output.resolve()
    overrides=json.loads(a.local_time_overrides.read_text()) if a.local_time_overrides else {}
    if out.exists():p.error('Use new output directory')
    refs={r['site_id']:r for r in csv.DictReader((ROOT/'data/reference/prc_spider_18_sites_verified.csv').open())}
    fig,axes=plt.subplots(1,len(a.sites),figsize=(6.5*len(a.sites),6),squeeze=False);results={}
    for ax,site in zip(axes.flat,a.sites):
        site_id,zone=SITES[site]
        folder=a.batch/site/'registration';report=json.loads((folder/'rotation_analysis.json').read_text())
        raw=next(r for r in report['inputs'] if Path(r['path']).name=='1.tif')
        explicit=overrides.get(site)
        if site=='Cottage' and a.cottage_local_time:explicit=a.cottage_local_time
        if site=='oxfordbussinesspark' and a.oxford_local_time:explicit=a.oxford_local_time
        stamp,origin=(explicit,'explicit user timestamp') if explicit else capture_time(Path(raw['path']))
        utc,offset=local_to_utc(stamp,zone)
        reference=refs[site_id];lat=float(reference['latitude']);lon=float(reference['longitude'])
        alt,az=lunar_track([utc],lat,lon);alt=float(alt[0]);az=float(az[0])
        with np.load(folder/'linear_stokes.npz',allow_pickle=False) as z:i=z['I']
        center=report['reference_center_xy_px'];r90=report['radius_60_px']*np.sqrt(2)
        y,x=np.indices(i.shape);radius=np.hypot(x-center[0],y-center[1]);sky=radius<r90
        view_theta=2*np.arcsin(np.clip(radius/(r90*np.sqrt(2)),0,1));view_alt=90-np.rad2deg(view_theta)
        flags=[];candidate=None;ring=None
        date=datetime.strptime(stamp,'%Y:%m:%d %H:%M:%S').strftime('%Y-%m-%d')
        if date!=reference['field_date']:flags.append('Image metadata date differs from reference field date; verify before interpretation')
        if 'fallback' in origin:flags.append('TIFF DateTime may be export timestamp')
        if alt<=0:flags.append('Moon predicted below horizon: timestamp/identity requires checking')
        else:
            ring=r90*np.sin(np.deg2rad((90-alt)/2))/np.sin(np.pi/4)
            permitted=sky&(abs(view_alt-alt)<=5)
            clean=np.nan_to_num(i)
            score=ndimage.gaussian_filter(clean,3)-ndimage.gaussian_filter(clean,20)
            score=np.where(permitted,score,-np.inf)
            iy,ix=np.unravel_index(np.argmax(score),score.shape);candidate=[float(ix),float(iy)]
            flags.append('Ring-constrained brightness candidate: Moon identity and image north not independently verified')
        ax.imshow(np.where(sky,i,np.nan),cmap='gray',vmin=np.nanpercentile(i[sky],5),vmax=np.nanpercentile(i[sky],99.9))
        if ring is not None:ax.add_patch(Circle(center,ring,fill=False,color='cyan',lw=1.2))
        if candidate is not None:ax.plot(*candidate,'rx',ms=9)
        ax.set_title(f'{site}\n{stamp} local → {utc} UTC\nMoon elevation {alt:.2f}deg; azimuth {az:.2f}deg',fontsize=9)
        ax.set_xlabel('Original pixel x');ax.set_ylabel('Original pixel y')
        results[site]=dict(site_id=site_id,local_time=stamp,timezone=zone,utc_offset_hours=offset,utc=utc,
            timestamp_source=origin,reference_field_date=reference['field_date'],moon_elevation_deg=alt,
            moon_compass_azimuth_deg=az,predicted_ring_radius_px=ring,candidate_pixel_xy=candidate,
            flags=flags,source_report_sha256=digest(folder/'rotation_analysis.json'))
        print(f'{site}: {utc} UTC, Moon elevation {alt:.2f}deg; '+ '; '.join(flags),flush=True)
    out.mkdir(parents=True);fig.suptitle('Cyan = ephemeris elevation ring; red = provisional candidate on ring\nCamera north unknown; inspect candidates and metadata dates before model interpretation')
    fig.tight_layout(rect=[0,0,1,.92]);fig.savefig(out/'moon_ephemeris_anchors.png',dpi=170);plt.close(fig)
    (out/'moon_ephemeris_anchors.json').write_text(json.dumps(dict(sites=results,
        status='EPHEMERIS_GEOMETRY_WITH_PROVISIONAL_IMAGE_BEARING',script_sha256=digest(__file__)),indent=2)+'\n')

if __name__=='__main__':main()
