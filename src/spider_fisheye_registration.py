"""Tilt-aware equisolid registration from image landmarks, not polarization."""
import numpy as np
from scipy import ndimage,optimize
from scipy.spatial.transform import Rotation


def rays(points,center,radius):
    delta=(np.asarray(points)-center)/radius
    r=np.linalg.norm(delta,axis=-1)
    theta=2*np.arcsin(np.clip(r/np.sqrt(2),0,.999999))
    factor=np.sin(theta)/np.maximum(r,1e-12)
    return np.concatenate([delta*factor[...,None],np.cos(theta)[...,None]],axis=-1)


def pixels(vectors,center,radius):
    vectors=np.asarray(vectors)
    rr=radius*np.sqrt(np.maximum(0,1-vectors[...,2]))
    norm=np.linalg.norm(vectors[...,:2],axis=-1)
    return center+vectors[...,:2]*(rr/np.maximum(norm,1e-12))[...,None]


def kabsch(source,target):
    u,_,vt=np.linalg.svd(source.T@target)
    sign=np.eye(3);sign[2,2]=np.linalg.det(vt.T@u.T)
    return vt.T@sign@u.T


def forward(points,mapping):
    if 'shared_tangential_p1_p2' in mapping:
        from spider_shared_fisheye_registration import forward as shared_forward
        return shared_forward(points,mapping)
    return pixels(rays(points,np.array(mapping['source_center_xy_px']),mapping['radius_px'])@
        np.array(mapping['ray_rotation_source_to_reference']).T,
        np.array(mapping['target_projection_center_xy_px']),mapping['radius_px'])


def fit_correspondences(source,target,center,radius):
    source=np.asarray(source);target=np.asarray(target);center=np.asarray(center)
    held=np.arange(len(source))%4==0
    s=source[~held];t=target[~held]
    if len(s)<8 or held.sum()<4:raise ValueError('Insufficient train/held-out correspondences')
    sr=rays(s,center,radius);tr=rays(t,center,radius)
    rng=np.random.default_rng(42);best_count=-1;initial=None
    for _ in range(300):
        ids=rng.choice(len(s),3,replace=False)
        rotation=kabsch(sr[ids],tr[ids])
        errors=np.linalg.norm(pixels(sr@rotation.T,center,radius)-t,axis=1)
        count=int((errors<2.5).sum())
        if count>best_count:best_count=count;initial=rotation
    start=np.r_[Rotation.from_matrix(initial).as_rotvec(),np.zeros(2)]
    def residual(q):
        rotation=Rotation.from_rotvec(q[:3]).as_matrix()
        errors=(pixels(sr@rotation.T,center+q[3:5],radius)-t).ravel()
        # Centre/tilt parameters can trade off. A stated weak pixel prior keeps
        # poorly identified centre offsets close to the inherited calibration.
        return np.r_[errors,.2*q[3:]]
    bounds=(np.r_[[-np.inf]*3,[-.05*radius]*2],np.r_[[np.inf]*3,[.05*radius]*2])
    result=optimize.least_squares(residual,start,bounds=bounds,loss='soft_l1',f_scale=1,max_nfev=200)
    rotation=Rotation.from_rotvec(result.x[:3]).as_matrix()
    mapping=dict(kind='equisolid_ray_rotation_with_bounded_projection_center_offsets',
        source_center_xy_px=center.tolist(),
        target_projection_center_xy_px=(center+result.x[3:5]).tolist(),radius_px=float(radius),
        ray_rotation_source_to_reference=rotation.tolist(),
        fitted_image_geometry_only=True,center_offset_bound_fraction=.05,center_offset_prior_weight=.2)
    error=np.linalg.norm(forward(source,mapping)-target,axis=1)
    inside=error[~held]<2.5;validation=error[held];confirmed=validation<3
    good=target[~held][inside]
    sectors=((np.arctan2(good[:,1]-center[1],good[:,0]-center[0])+np.pi)*8/(2*np.pi)).astype(int)%8
    span=float(np.linalg.norm(np.ptp(good,axis=0))) if len(good) else 0.
    boundary=bool(np.any(abs(result.x[3:])>.0495*radius))
    quality=dict(training_matches=int((~held).sum()),training_inliers=int(inside.sum()),
        training_median_px=float(np.median(error[~held][inside])) if inside.any() else None,
        heldout_matches=int(held.sum()),heldout_inliers=int(confirmed.sum()),
        heldout_inlier_fraction=float(confirmed.mean()),heldout_median_all_matches_px=float(np.median(validation)),
        heldout_inlier_p95_px=float(np.percentile(validation[confirmed],95)) if confirmed.any() else None,
        occupied_reference_sectors=int(len(np.unique(sectors))),landmark_span_px=span,
        optimizer_success=bool(result.success),center_offset_limit_reached=boundary,
        projection_center_offsets_xy_px=result.x[3:].tolist())
    quality['numerical_pass']=bool(result.success and not boundary and inside.sum()>=8 and confirmed.sum()>=4
        and confirmed.mean()>=.6 and np.median(validation)<2 and len(np.unique(sectors))>=3 and span>.5*radius)
    return mapping,quality,error


def warp(image,mapping,shape=None):
    """Inverse mapping; interpolate raw exposure once, unsupported pixels NaN."""
    if 'shared_tangential_p1_p2' in mapping:
        from spider_shared_fisheye_registration import warp as shared_warp
        return shared_warp(image,mapping,shape)
    shape=image.shape if shape is None else shape
    y,x=np.indices(shape,dtype=float)
    target=np.stack([x,y],axis=-1)
    source_vectors=(rays(target,np.array(mapping['target_projection_center_xy_px']),mapping['radius_px'])@
        np.array(mapping['ray_rotation_source_to_reference']))
    src=pixels(source_vectors,np.array(mapping['source_center_xy_px']),mapping['radius_px'])
    result=ndimage.map_coordinates(np.asarray(image,dtype=float),[src[...,1],src[...,0]],order=1,
        mode='constant',cval=np.nan,prefilter=False)
    supported=np.linalg.norm(target-np.array(mapping['target_projection_center_xy_px']),axis=-1)<1.03*mapping['radius_px']
    return np.where(supported,result,np.nan)

