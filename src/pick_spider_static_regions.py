"""Select stationary reference landmark regions for geometric registration."""
import argparse
import base64
import hashlib
import io
import json
from pathlib import Path
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
from matplotlib.image import imsave


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


HTML = r"""<!doctype html><html lang="en"><meta charset="utf-8"><title>Static landmark regions</title>
<style>body{font:16px system-ui;margin:20px;background:#eee}button,input{font:inherit;margin:5px}#view{height:75vh;overflow:auto;background:#111}canvas{display:block;cursor:crosshair}</style>
<h2 id="heading"></h2><p>Outline narrow regions around stationary skyline, trees, buildings or other ground landmarks.
Include skyline edges, but avoid open star fields, the Moon, light trails and moving vegetation where possible.
Select separated regions around at least three sides of the image. These regions constrain registration only, not the DoLP analysis field.</p>
<p>Click polygon vertices, then Finish region. Repeat for additional regions. Download JSON when done.</p>
<button id="finish">Finish region</button><button id="undo">Undo vertex</button><button id="remove">Remove last region</button><button id="save">Download regions JSON</button>
<label>Zoom <input id="zoom" type="range" min="25" max="250" value="75"></label>
<label>Brightness <input id="brightness" type="range" min="50" max="600" value="150"></label>
<p id="status"></p><div id="view"><canvas id="canvas"></canvas></div>
<script>
const metadata=__METADATA__,image=new Image(),canvas=document.getElementById('canvas'),ctx=canvas.getContext('2d');
let current=[],regions=[];
document.getElementById('heading').textContent=metadata.site+' — select stationary landmark regions';
function draw(){ctx.clearRect(0,0,canvas.width,canvas.height);ctx.drawImage(image,0,0);
 for(const [points,closed] of [...regions.map(r=>[r,true]),[current,false]]){
  if(!points.length)continue;ctx.beginPath();ctx.moveTo(...points[0]);for(const p of points.slice(1))ctx.lineTo(...p);
  if(closed){ctx.closePath();ctx.fillStyle='rgba(0,255,90,.18)';ctx.fill();}
  ctx.strokeStyle=closed?'lime':'yellow';ctx.lineWidth=3;ctx.stroke();
  for(const p of points){ctx.beginPath();ctx.arc(...p,4,0,2*Math.PI);ctx.fillStyle='yellow';ctx.fill();}
 }
 document.getElementById('status').textContent=regions.length+' completed regions; '+current.length+' current vertices';
}
function zoom(){const f=+document.getElementById('zoom').value/100;canvas.style.width=image.width*f+'px';canvas.style.height=image.height*f+'px';}
image.onload=()=>{canvas.width=image.width;canvas.height=image.height;zoom();draw();};
image.src='data:image/png;base64,__IMAGE__';
canvas.onclick=e=>{const r=canvas.getBoundingClientRect();current.push([(e.clientX-r.left)*canvas.width/r.width,(e.clientY-r.top)*canvas.height/r.height]);draw();};
document.getElementById('finish').onclick=()=>{if(current.length<3){alert('At least three vertices required.');return;}regions.push(current);current=[];draw();};
document.getElementById('undo').onclick=()=>{current.pop();draw();};
document.getElementById('remove').onclick=()=>{regions.pop();draw();};
document.getElementById('zoom').oninput=zoom;
document.getElementById('brightness').oninput=e=>canvas.style.filter='brightness('+e.target.value/100+')';
document.getElementById('save').onclick=()=>{
 if(current.length){alert('Finish or undo the current polygon first.');return;}
 if(!regions.length){alert('Select at least one region.');return;}
 const data={...metadata,regions_xy:regions,purpose:'Investigator-selected stationary reference landmark regions; celestial sources excluded where possible'};
 const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'})),a=document.createElement('a');
 a.href=url;a.download=metadata.site+'_static_landmark_regions.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),2000);
};
</script></html>"""



def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-lock',type=Path,required=True)
    p.add_argument('--site',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    lock=json.loads(a.reference_lock.read_text())['sites'][a.site]
    candidates=[r for r in lock['raw_inputs'] if Path(r['path']).stem=='1']
    if len(candidates)!=1:p.error('Exactly one reference exposure 1 required')
    record=candidates[0];path=Path(record['path'])
    if digest(path)!=record['sha256']:p.error('Reference TIFF checksum changed')
    image=tifffile.imread(path)
    if image.ndim!=2:p.error('Requires prepared single-channel reference TIFF')
    finite=image[np.isfinite(image)]
    if not finite.size:p.error('No finite reference pixels')
    lo,hi=np.percentile(finite,[5,99.9])
    if hi<=lo:p.error('Reference has no display contrast')
    display=np.clip(np.nan_to_num((image-lo)/(hi-lo)),0,1)
    buffer=io.BytesIO()
    imsave(buffer,display,cmap='gray',vmin=0,vmax=1,format='png',origin='upper')
    metadata=dict(site=a.site,reference_path=str(path),reference_sha256=digest(path),
        source_lock_path=str(a.reference_lock),source_lock_sha256=digest(a.reference_lock),
        image_shape_yx=list(image.shape),reference_center_xy_px=lock['reference_center_xy_px'],
        rim_radius_px=lock['rim_radius_px'],coordinate_convention='Original reference pixels; x right, y down; no rotation or resampling',
        status='Picked point candidates only; no celestial orientation inferred',
        picker_script_sha256=digest(__file__))
    # Escape HTML-sensitive characters in JSON even if a source path contains them.
    serialized=json.dumps(metadata).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    html=HTML.replace('__METADATA__',serialized).replace('__IMAGE__',base64.b64encode(buffer.getvalue()).decode())
    a.output.mkdir(parents=True)
    (a.output/'pick_static_regions.html').write_text(html)
    (a.output/'reference_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Open '+str(a.output/'pick_static_regions.html'))
    print('Download regions JSON before closing the browser.')


if __name__=='__main__':main()
