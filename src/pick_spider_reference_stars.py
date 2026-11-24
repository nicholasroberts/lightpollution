"""Create an offline browser picker for stars/landmarks in a locked reference TIFF."""
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


HTML = r"""<!doctype html>
<html lang="en"><meta charset="utf-8"><title>Reference star picker</title>
<style>
body{font:16px system-ui;margin:20px;background:#eee;color:#222}
button,input{font:inherit;margin:4px} #viewport{height:72vh;overflow:auto;background:#111}
canvas{display:block;cursor:crosshair;transform-origin:top left}
table{border-collapse:collapse} td,th{padding:4px 12px;text-align:left}
</style>
<h2 id="heading"></h2>
<p>Click the centres of stars or other identifiable points. Pick at least two stars; additional stars help resolve identities.
Coordinates are in the original reference image, x right / y down. No star identity is assumed.</p>
<label>Next label <input id="label" placeholder="P1, P2 … (or a known name)"></label>
<button id="undo">Undo last</button>
<button id="save">Download coordinates JSON</button>
<button id="preview">Download marked PNG</button><br>
<label>Zoom <input id="zoom" type="range" min="25" max="250" value="75"> <span id="ztext">75%</span></label>
<label>Brightness <input id="brightness" type="range" min="50" max="400" value="100"></label>
<p>Use zoom and the scroll bars for accurate clicks. Download the JSON before closing; selections are not automatically saved.</p>
<div id="viewport"><canvas id="canvas"></canvas></div>
<table><thead><tr><th>Point</th><th>x</th><th>y</th></tr></thead><tbody id="points"></tbody></table>
<script>
const metadata=__METADATA__;
const image=new Image();
const canvas=document.getElementById('canvas'),ctx=canvas.getContext('2d');
const points=[],table=document.getElementById('points');
document.getElementById('heading').textContent=metadata.site+' — reference exposure point picker';
function redraw(){
 ctx.clearRect(0,0,canvas.width,canvas.height);
 ctx.drawImage(image,0,0);
 ctx.strokeStyle='#ff4545';ctx.fillStyle='#ff4545';ctx.lineWidth=2;ctx.font='20px sans-serif';
 for(const p of points){
  const [x,y]=p.pixel_xy;
  ctx.beginPath();ctx.arc(x,y,9,0,2*Math.PI);ctx.stroke();
  ctx.beginPath();ctx.moveTo(x-14,y);ctx.lineTo(x+14,y);ctx.moveTo(x,y-14);ctx.lineTo(x,y+14);ctx.stroke();
  ctx.fillText(p.label,x+16,y-12);
 }
 table.replaceChildren();
 for(const p of points){
  const row=document.createElement('tr');
  for(const v of [p.label,...p.pixel_xy.map(v=>v.toFixed(2))]){
   const cell=document.createElement('td');cell.textContent=v;row.appendChild(cell);
  } table.appendChild(row);
 }
}
function zoom(){
 const factor=Number(document.getElementById('zoom').value)/100;
 canvas.style.width=(image.width*factor)+'px';canvas.style.height=(image.height*factor)+'px';
 document.getElementById('ztext').textContent=Math.round(factor*100)+'%';
}
image.onload=()=>{canvas.width=image.width;canvas.height=image.height;zoom();redraw();};
image.src='data:image/png;base64,__IMAGE__';
canvas.addEventListener('click',event=>{
 const rect=canvas.getBoundingClientRect();
 const x=(event.clientX-rect.left)*canvas.width/rect.width;
 const y=(event.clientY-rect.top)*canvas.height/rect.height;
 const label=document.getElementById('label').value.trim()||'P'+(points.length+1);
 points.push({label,pixel_xy:[x,y],identity_status:'user-selected point; astronomical identity unverified'});
 document.getElementById('label').value='';redraw();
});
document.getElementById('undo').onclick=()=>{points.pop();redraw();};
document.getElementById('zoom').oninput=zoom;
document.getElementById('brightness').oninput=event=>{
 canvas.style.filter='brightness('+Number(event.target.value)/100+')';
};
function download(blob,name){
 const url=URL.createObjectURL(blob),link=document.createElement('a');
 link.href=url;link.download=name;link.click();setTimeout(()=>URL.revokeObjectURL(url),2000);
}
document.getElementById('save').onclick=()=>{
 if(points.length<2){alert('Select at least two points before saving.');return;}
 download(new Blob([JSON.stringify({...metadata,points},null,2)],{type:'application/json'}),
 metadata.site+'_reference_star_points.json');
};
document.getElementById('preview').onclick=()=>canvas.toBlob(blob=>{
 download(blob,metadata.site+'_reference_star_points.png');
},'image/png');
</script></html>
"""


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
    (a.output/'pick_stars.html').write_text(html)
    (a.output/'reference_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Open '+str(a.output/'pick_stars.html'))
    print('Download coordinates JSON and marked PNG after selecting points.')


if __name__=='__main__':main()
