import argparse,json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter,minimum_filter,map_coordinates
from scipy.optimize import root

R=655.544999673502
def ray(xy,c):
 v=np.asarray(xy)-c;r=np.linalg.norm(v,axis=-1);t=2*np.arcsin(r/(R*np.sqrt(2)));a=np.arctan2(-v[...,1],v[...,0]);return np.stack([np.sin(t)*np.cos(a),np.sin(t)*np.sin(a),np.cos(t)],-1)
def angle(a,b):return float(np.degrees(np.arccos(np.clip(np.dot(a,b),-1,1))))
parser=argparse.ArgumentParser(description='Measure neutral points in the two reviewed uploaded north maps; fixed diagnostic windows, no model fit')
parser.add_argument('--maps',type=Path,required=True,help='Contains Cottage north_stokes.npz and museum north_stokes(1).npz')
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
if args.output.exists():parser.error('Use fresh output JSON')
out=[]
for file,site,el,az,boxes in [('north_stokes.npz','Cottage',12.31370836118817,191.15960520703004,[(850,1200,800,1080),(700,1120,90,300)]),('north_stokes(1).npz','museum',41.19831511184027,187.1466130743674,[(1150,1430,250,530),(650,900,640,900)])]:
 z=np.load(args.maps/file);c=z['center_xy'];i=z['I'];q=z['Q'];u=z['U'];valid=z['valid'] & np.isfinite(i)&np.isfinite(q)&np.isfinite(u)
 # Smooth Stokes in small diagnostic windows to suppress stars/pixel noise; never smooth angles.
 for sigma in [4,8,12]:
  w=gaussian_filter(valid.astype(float),sigma)
  ii,qq,uu=[gaussian_filter(np.where(valid,v,0),sigma)/np.maximum(w,1e-15) for v in (i,q,u)]
  d=np.hypot(qq,uu)/ii
  moon=np.array([-np.cos(np.radians(el))*np.sin(np.radians(az)),np.cos(np.radians(el))*np.cos(np.radians(az)),np.sin(np.radians(el))])
  coords=[]
  for box in boxes:
   x0,x1,y0,y1=box;patch=d[y0:y1,x0:x1];idx=np.unravel_index(np.nanargmin(patch),patch.shape);xy=np.array([idx[1]+x0,idx[0]+y0],float)
   # Root of smoothed Q,U local field, only accepted inside chosen diagnostic window.
   def fun(p):return [float(map_coordinates(v,[[p[1]],[p[0]]],order=1)[0]) for v in (qq,uu)]
   sol=root(fun,xy);p=sol.x if sol.success and x0<sol.x[0]<x1 and y0<sol.x[1]<y1 else xy
   rr=ray(p,c);alt=float(np.degrees(np.arcsin(rr[2])));azi=float(np.degrees(np.arctan2(-rr[0],rr[1]))%360)
   coords.append(dict(xy=p.tolist(),root=bool(sol.success and np.linalg.norm(fun(p))<1e-4),elevation=alt,azimuth=azi,moon_distance=angle(rr,moon),dolp=float(map_coordinates(d,[[p[1]],[p[0]]],order=1)[0])))
  sep=angle(ray(coords[0]['xy'],c),ray(coords[1]['xy'],c))
  print(site,sigma,coords,'visible-points separation',sep)
  out.append(dict(site=site,sigma_px=sigma,points=coords,visible_points_separation_deg=sep))
args.output.write_text(json.dumps(out,indent=2))
