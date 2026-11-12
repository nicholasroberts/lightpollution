"""Provisional fisheye geometry audit; never applies calibration or changes TIFFs.

Run on original Cottage TIFFs. Fits the outer illuminated boundary at three
thresholds, independently for every exposure. Trees, vignetting, tilt and crop
can bias the result: geometry_overlay.png MUST be reviewed. An optical rim is
not a measured horizon. Nominal equisolid conversion assumes rim = theta 90deg.
60deg means radial angle, not a 60deg-diameter field.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy import ndimage, optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import tifffile


def fit_boundary(image, fraction):
    image = np.asarray(image, dtype=float)
    if image.ndim != 2 or not np.isfinite(image).all():
        raise ValueError('Expected a finite 2D channel')
    h, w = image.shape
    baseline = np.percentile(image, 5)
    scale = np.percentile(image, 95) - baseline
    if scale <= 0:
        raise ValueError('Insufficient intensity range')
    mask = image > baseline + fraction * scale
    labels, n = ndimage.label(mask)
    if not n:
        raise ValueError('No illuminated component')
    sizes = np.bincount(labels.ravel()); sizes[0] = 0
    mask = labels == sizes.argmax()
    mask = ndimage.binary_fill_holes(mask)
    edge = mask & ~ndimage.binary_erosion(mask)
    y, x = np.nonzero(edge)
    # Select the outermost boundary per angular sector to reduce star/hole bias.
    cx, cy = (w-1)/2, (h-1)/2
    angle = np.arctan2(y-cy, x-cx)
    sector = np.floor((angle + np.pi)/(2*np.pi)*360).astype(int).clip(0,359)
    radius = np.hypot(x-cx, y-cy)
    indices = [ix[np.argmax(radius[ix])] for k in range(360)
               if len(ix := np.flatnonzero(sector == k))]
    x, y = x[indices], y[indices]
    if len(x) < 180:
        raise ValueError('Insufficient boundary coverage')
    def residual(p):
        return np.hypot(x-p[0], y-p[1])-p[2]
    fit = optimize.least_squares(residual, [cx,cy,np.median(radius)],
        loss='soft_l1', f_scale=3,
        bounds=([0,0,min(h,w)*.2],[w-1,h-1,min(h,w)*.65]))
    cx,cy,r = fit.x
    err = residual(fit.x)
    return dict(center_x_px=float(cx),center_y_px=float(cy),rim_radius_px=float(r),
        r60_equisolid_nominal_px=float(r/np.sqrt(2)),
        threshold_fraction=fraction,median_abs_residual_px=float(np.median(abs(err))),
        p90_abs_residual_px=float(np.percentile(abs(err),90)),
        sectors=len(x),circle_within_frame=bool(cx-r>0 and cy-r>0 and cx+r<w-1 and cy+r<h-1)), x,y


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--channel',type=int,default=2,choices=[0,1,2])
    a=p.parse_args()
    files=sorted([f for f in a.folder.iterdir() if f.suffix.lower() in ('.tif','.tiff') and f.stem.isdigit()], key=lambda f:int(f.stem))
    if not files: p.error('No numbered TIFFs found')
    a.output.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(int(np.ceil(len(files)/4)),4,figsize=(16,4*np.ceil(len(files)/4)),squeeze=False)
    records=[]
    for f,ax in zip(files,axes.flat):
        raw=tifffile.imread(f)
        if raw.ndim==3 and raw.shape[-1]==3: image=raw[...,a.channel]
        elif raw.ndim==2: image=raw
        else: raise ValueError(f'Unsupported shape {raw.shape}: {f}')
        ax.imshow(image,cmap='gray',vmin=np.percentile(image,1),vmax=np.percentile(image,99.5))
        record=dict(file=str(f.resolve()),shape=list(raw.shape),dtype=str(raw.dtype),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),fits=[])
        for threshold,color in zip([.02,.05,.10],['cyan','yellow','magenta']):
            try:
                fit,x,y=fit_boundary(image,threshold)
                record['fits'].append(fit)
                c=(fit['center_x_px'],fit['center_y_px'])
                ax.add_patch(Circle(c,fit['rim_radius_px'],fill=False,color=color,lw=.8))
                ax.add_patch(Circle(c,fit['r60_equisolid_nominal_px'],fill=False,color=color,lw=.8,ls='--'))
                ax.plot(*c,'+',color=color,ms=5)
                ax.scatter(x[::6],y[::6],s=1,color=color)
            except ValueError as e: record['fits'].append(dict(threshold_fraction=threshold,error=str(e)))
        ax.set_title(f.name); ax.axis('off'); records.append(record)
    for ax in list(axes.flat)[len(files):]: ax.axis('off')
    fig.suptitle('PROVISIONAL: solid = fitted rim; dashed = nominal 60deg radial field\ncyan/yellow/magenta = 2/5/10% thresholds; review boundary, centre and crop')
    fig.tight_layout(rect=[0,0,1,.94]); fig.savefig(a.output/'geometry_overlay.png',dpi=160); plt.close(fig)
    report=dict(status='provisional_requires_visual_review_not_applied',
        assumptions=['Nominal equisolid projection r(theta)=2f sin(theta/2)',
        'Fitted illuminated boundary treated provisionally as 90deg optical rim',
        'Camera optical axis assumed zenith for zenith-centred interpretation',
        'No angular response weighting; no registration or Stokes calculation performed'],
        cautions=['Threshold spread is a diagnostic, not a confidence interval',
        'Dark ground, vignetting and obstructions can move the detected boundary',
        'Do not transfer to other crops, sizes or geometries without checking'],
        channel=a.channel,records=records,raw_images_modified=False)
    (a.output/'geometry_calibration.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'Wrote {a.output}/geometry_overlay.png and geometry_calibration.json; review before use')

if __name__=='__main__': main()
