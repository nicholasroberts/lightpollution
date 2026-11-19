"""Create exact-pixel Moon selection and saved-landmark residual diagnostics.

No exposure transforms or astronomical coordinates are changed by this tool.
"""
import argparse
import hashlib
import html
import json
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_page(a):
    from PIL import Image
    import tifffile
    rows=json.loads((a.registration/'registration_batch_summary.json').read_text())['sites']
    indexed={r['site']:r for r in rows}
    cards=[];records={}
    # Validate sources before writing the fresh review directory.
    prepared=[]
    for site in a.sites:
        if site not in indexed:raise ValueError('Unknown site: '+site)
        path=a.raw_root/site/'1.tif'
        manifest_path=path.parent/'linear_raw_manifest.json'
        manifest=json.loads(manifest_path.read_text())
        frame=next(f for f in manifest['frames'] if f['exposure']==1)
        if digest(path)!=frame['output_sha256']:raise ValueError('Reference TIFF hash mismatch: '+site)
        row=indexed[site]
        if row.get('lock_path') and digest(a.registration/row['lock_path'])!=row['lock_sha256']:
            raise ValueError('Registration lock hash mismatch: '+site)
        pixels=tifffile.imread(path)
        if pixels.ndim!=2 or not np.isfinite(pixels).all():raise ValueError('Expected finite single-channel reference')
        prepared.append((site,path,pixels,manifest_path,row))
    a.output.mkdir(parents=True,exist_ok=False)
    for site,path,pixels,manifest_path,row in prepared:
        positive=pixels[pixels>0]
        if not positive.size:raise ValueError('No positive reference pixels: '+site)
        lo,hi=np.percentile(positive,[1,99.9])
        display=np.log1p(100*np.clip((pixels-lo)/max(hi-lo,1e-12),0,1))/np.log(101)
        Image.fromarray(np.uint8(255*display)).save(a.output/(site+'.png'))
        records[site]=dict(reference_tiff_sha256=digest(path),linear_raw_manifest_sha256=digest(manifest_path),
            image_shape_yx=list(pixels.shape),reference_exposure='1.tif',
            registration_lock_sha256=row.get('lock_sha256'))
        cards.append(f'<section data-site="{html.escape(site)}"><h2>{html.escape(site)}</h2><p class="position">No selection</p><button class="clear">Clear selection</button><label>Note <input value="Investigator identified Moon in linear green reference 1.tif"></label><div class="view"><img src="{site}.png" draggable="false"><span class="marker">+</span></div></section>')
    script=r'''
const references=REFERENCES;
const points={};
for(const card of document.querySelectorAll('section')){
 const view=card.querySelector('.view'),img=card.querySelector('img'),mark=card.querySelector('.marker');
 img.onclick=e=>{const b=img.getBoundingClientRect();const x=Math.max(0,Math.min(img.naturalWidth-1,Math.round((e.clientX-b.left)/b.width*img.naturalWidth-.5)));const y=Math.max(0,Math.min(img.naturalHeight-1,Math.round((e.clientY-b.top)/b.height*img.naturalHeight-.5)));points[card.dataset.site]=[x,y];mark.style.left=((x+.5)/img.naturalWidth*100)+'%';mark.style.top=((y+.5)/img.naturalHeight*100)+'%';mark.style.display='block';card.querySelector('.position').textContent='Moon pixel x='+x+', y='+y;};
 card.querySelector('.clear').onclick=()=>{delete points[card.dataset.site];mark.style.display='none';card.querySelector('.position').textContent='No selection';};
}
document.querySelector('#zoom').onchange=e=>document.body.classList.toggle('native',e.target.checked);
document.querySelector('#save').onclick=()=>{const result={};for(const card of document.querySelectorAll('section')){const s=card.dataset.site;if(points[s]){const note=card.querySelector('input').value.trim();if(!note){alert('Supply an identification note for '+s);return;}result[s]={moon_pixel_xy:points[s],note,...references[s]};}}if(!Object.keys(result).length){alert('Select at least one Moon');return;}const url=URL.createObjectURL(new Blob([JSON.stringify(result,null,2)],{type:'application/json'}));const link=document.createElement('a');link.href=url;link.download='moon_pixel_overrides.json';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
'''.replace('REFERENCES',json.dumps(records))
    page='''<!doctype html><meta charset="utf-8"><title>Identify reference Moon pixels</title><style>body{font:17px system-ui;margin:24px}header{position:sticky;top:0;background:white;padding:15px;z-index:2}section{border-top:1px solid #aaa;padding:15px 0}.view{position:relative;width:100%;max-width:2000px}.view img{width:100%;display:block;cursor:crosshair}.marker{display:none;position:absolute;transform:translate(-50%,-50%);color:red;font:bold 28px monospace;pointer-events:none}.native .view{width:2000px;max-width:none}input{width:550px}header input{width:auto}</style><header><button id="save">Download Moon coordinates</button><label><input id="zoom" type="checkbox"> Enlarge for precise selection</label></header><h1>Identify the Moon in reference 1.tif</h1><p>Click the Moon centre. These are full-canvas previews with the original pixel dimensions and no rotation, crop or reflection. The contrast stretch is for display only. Blank sites are omitted from export. Download before closing; selections are not saved automatically. No transforms or sky models are fitted here.</p>'''+''.join(cards)+'<script>'+script+'</script>'
    (a.output/'moon_pixels.html').write_text(page,encoding='utf-8')
    return indexed


def residual_diagnostics(a,indexed):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    result={}
    for site in a.diagnose:
        row=indexed[site]
        if not row.get('lock_path'):continue
        path=a.registration/row['lock_path']
        if digest(path)!=row['lock_sha256']:raise ValueError('Changed registration lock: '+site)
        lock=json.loads(path.read_text())['sites'][site]
        fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
        failed=[]
        for t in lock.get('transforms',[]):
            q=t['quality'];name=t['filename'];p=t.get('landmark_correspondences')
            if not p:continue
            held=np.asarray(p['heldout'],dtype=bool);err=np.asarray(p['residual_px'])
            vals=err[held]
            axes[0].plot(np.sort(vals),np.linspace(0,1,len(vals)),label=name)
            if not q.get('numerical_pass'):
                xy=np.asarray(p['reference_xy_px'])[held]
                im=axes[1].scatter(xy[:,0],xy[:,1],c=vals,vmin=0,vmax=5,s=12)
                fig.colorbar(im,ax=axes[1],label='Held-out residual (px)')
                failed.append(dict(filename=name,quality=q))
        axes[0].axvline(2,color='black',ls='--');axes[0].set(xlim=(0,6),xlabel='Held-out residual (px)',ylabel='Cumulative fraction');axes[0].legend(fontsize=7)
        axes[1].invert_yaxis();axes[1].set_aspect('equal');axes[1].set_title('Failed exposures: held-out landmark locations')
        fig.suptitle(site+' — saved registration residuals; no refit or relaxed thresholds')
        fig.savefig(a.output/(site+'_residuals.png'),dpi=160);plt.close(fig)
        result[site]=dict(lock_sha256=digest(path),failed_exposures=failed)
    (a.output/'registration_residual_diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registration',type=Path,required=True)
    p.add_argument('--raw-root',type=Path,required=True,help='Prepared linear green TIFF batch root')
    p.add_argument('--sites',nargs='+',default=['WClawn','oxfordbussinesspark','ChurchHill','PennCC','Waitrose','minifactory'])
    p.add_argument('--diagnose',nargs='*',default=['minifactory'])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try:
        indexed=make_page(a);residual_diagnostics(a,indexed)
    except (OSError,ValueError,KeyError,StopIteration) as exc:p.error(str(exc))
    print('Open '+str(a.output/'moon_pixels.html'))


if __name__=='__main__':main()
