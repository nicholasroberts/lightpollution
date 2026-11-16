"""Shared image-geometry correction; never fits Stokes, Moon or model values.

A bounded two-coefficient tangential correction is common to every exposure.
Only SO(3) attitude and a bounded projection-centre offset vary per exposure.
This is an empirical registration refinement, not independent lens calibration.
Whole 96-pixel reference cells are held out before any geometry fitting.
"""
import numpy as np
from scipy import ndimage, optimize
from scipy.spatial.transform import Rotation
from spider_fisheye_registration import rays, pixels, fit_correspondences


def heldout_cells(target):
    cells=np.floor(np.asarray(target)/96).astype(int)
    return (cells[:,0]+2*cells[:,1])%4==0


def undistort(points,center,radius,coefficients):
    d=(np.asarray(points)-center)/radius
    x,y=d[...,0],d[...,1];rr=x*x+y*y;p1,p2=coefficients
    correction=np.stack([2*p1*x*y+p2*(rr+2*x*x),p1*(rr+2*y*y)+2*p2*x*y],axis=-1)
    return center+radius*(d+correction)


def distort(points,center,radius,coefficients):
    result=np.array(points,dtype=float,copy=True)
    for _ in range(12):
        result+=points-undistort(result,center,radius,coefficients)
    return result


def predict(source,pose,center,radius,coefficients):
    u=undistort(source,center,radius,coefficients)
    v=pixels(rays(u,center,radius)@Rotation.from_rotvec(pose[:3]).as_matrix().T,center+pose[3:5],radius)
    return distort(v,center,radius,coefficients)


def inverse_points(target,mapping):
    center=np.asarray(mapping['source_center_xy_px']);radius=mapping['radius_px']
    coefficients=mapping['shared_tangential_p1_p2'];tc=np.asarray(mapping['target_projection_center_xy_px'])
    u=undistort(target,center,radius,coefficients)
    v=pixels(rays(u,tc,radius)@np.asarray(mapping['ray_rotation_source_to_reference']),center,radius)
    return distort(v,center,radius,coefficients)


def forward(source,mapping):
    center=np.asarray(mapping['source_center_xy_px']);radius=mapping['radius_px']
    pose=np.r_[Rotation.from_matrix(mapping['ray_rotation_source_to_reference']).as_rotvec(),
               np.asarray(mapping['target_projection_center_xy_px'])-center]
    return predict(source,pose,center,radius,mapping['shared_tangential_p1_p2'])


def warp(image,mapping,shape=None):
    shape=image.shape if shape is None else shape
    y,x=np.indices(shape,dtype=float);target=np.stack([x,y],axis=-1)
    center=np.asarray(mapping['source_center_xy_px']);radius=mapping['radius_px']
    supported=np.linalg.norm(target-center,axis=-1)<1.03*radius
    # Keep the nonlinear inverse within its checked circular domain.
    source=inverse_points(target[supported],mapping)
    values=ndimage.map_coordinates(np.asarray(image,dtype=float),[source[:,1],source[:,0]],
        order=1,mode='constant',cval=np.nan,prefilter=False)
    result=np.full(shape,np.nan);result[supported]=values
    return result


def fit_site(correspondences,center,radius,calibration_exposures=None):
    """Return per-frame (mapping, quality, residuals) and shared diagnostics.

    For 12 exposures only 4 and 7 determine shared coefficients; 10--12
    therefore check transfer to an analyzer group not used for that fit.
    The same reference-cell holdout is excluded from both shared and pose fits.
    """
    center=np.asarray(center,dtype=float);radius=float(radius)
    if calibration_exposures is None:
        calibration_exposures=next((pair for pair in ([4,7],[4,10],[7,10],[2,3])
            if all(n in correspondences for n in pair)),None)
        if calibration_exposures is None:raise ValueError('Two distinct analyzer groups are required for shared geometry')
    if any(n not in correspondences for n in calibration_exposures):
        raise ValueError('Missing an explicit shared-geometry training exposure')
    records={};offset_bound=.05*radius
    for n,(source,target) in correspondences.items():
        source=np.asarray(source);target=np.asarray(target);held=heldout_cells(target)
        if (~held).sum()<16 or held.sum()<8:
            raise ValueError('Insufficient spatially separate train/holdout landmarks')
        # Initializer sees training coordinates only.
        initial,_,_=fit_correspondences(source[~held],target[~held],center,radius)
        pose=np.r_[Rotation.from_matrix(initial['ray_rotation_source_to_reference']).as_rotvec(),
                   np.asarray(initial['target_projection_center_xy_px'])-center]
        records[n]=(source,target,held,pose)
    start=np.r_[np.zeros(2),np.concatenate([records[n][3] for n in calibration_exposures])]
    def residual(v):
        chunks=[]
        for j,n in enumerate(calibration_exposures):
            source,target,held,_=records[n];pose=v[2+5*j:7+5*j]
            chunks.extend([(predict(source[~held],pose,center,radius,v[:2])-target[~held]).ravel(),.2*pose[3:]])
        return np.r_[np.concatenate(chunks),10*v[:2]]
    count=len(calibration_exposures)
    lower=np.r_[[-.01]*2,np.tile(np.r_[[-np.inf]*3,[-offset_bound]*2],count)]
    upper=np.r_[[.01]*2,np.tile(np.r_[[np.inf]*3,[offset_bound]*2],count)]
    joint=optimize.least_squares(residual,start,bounds=(lower,upper),loss='soft_l1',f_scale=1,max_nfev=300)
    coefficients=joint.x[:2]
    # Check inverse and local handedness over the full supported circular field.
    grid=np.linspace(-1.03,1.03,41);gx,gy=np.meshgrid(grid,grid);inside=gx*gx+gy*gy<=1.03**2
    gx,gy=gx[inside],gy[inside];p1,p2=coefficients
    determinants=(1+2*p1*gy+6*p2*gx)*(1+6*p1*gy+2*p2*gx)-(2*p1*gx+2*p2*gy)**2
    points=center+radius*np.c_[gx,gy]
    inverse_error=float(np.max(np.linalg.norm(distort(undistort(points,center,radius,coefficients),center,radius,coefficients)-points,axis=1)))
    shared_ok=bool(joint.success and np.all(abs(coefficients)<.0099) and determinants.min()>.5 and inverse_error<1e-5)
    shared=dict(kind='shared_tangential_image_geometry',coefficients_p1_p2=coefficients.tolist(),
        calibration_exposures=list(calibration_exposures),coefficient_bound=.01,
        reference_cell_size_px=96,heldout_rule='(floor(x/96)+2*floor(y/96)) % 4 == 0',
        optimizer_success=bool(joint.success),minimum_correction_jacobian_determinant=float(determinants.min()),
        inverse_roundtrip_max_px=inverse_error,numerical_geometry_check=shared_ok,
        independent_lens_calibration=False,polarization_or_astronomy_used_for_fit=False)
    fitted={}
    for n,(source,target,held,start_pose) in records.items():
        def pose_residual(p):
            return np.r_[(predict(source[~held],p,center,radius,coefficients)-target[~held]).ravel(),.2*p[3:]]
        fit=optimize.least_squares(pose_residual,start_pose,
            bounds=(np.r_[[-np.inf]*3,[-offset_bound]*2],np.r_[[np.inf]*3,[offset_bound]*2]),
            loss='soft_l1',f_scale=1,max_nfev=200)
        errors=np.linalg.norm(predict(source,fit.x,center,radius,coefficients)-target,axis=1)
        training_inliers=errors[~held]<2.5;validation=errors[held];confirmed=validation<3
        good=target[~held][training_inliers]
        sectors=((np.arctan2(good[:,1]-center[1],good[:,0]-center[0])+np.pi)*8/(2*np.pi)).astype(int)%8
        span=float(np.linalg.norm(np.ptp(good,axis=0))) if len(good) else 0.
        bound=bool(np.any(abs(fit.x[3:])>.0495*radius))
        quality=dict(training_matches=int((~held).sum()),training_inliers=int(training_inliers.sum()),
            training_median_px=float(np.median(errors[~held][training_inliers])) if training_inliers.any() else None,
            heldout_matches=int(held.sum()),heldout_inliers=int(confirmed.sum()),heldout_inlier_fraction=float(confirmed.mean()),
            heldout_median_all_matches_px=float(np.median(validation)),
            heldout_inlier_p95_px=float(np.percentile(validation[confirmed],95)) if confirmed.any() else None,
            occupied_reference_sectors=len(np.unique(sectors)),landmark_span_px=span,
            optimizer_success=bool(fit.success),center_offset_limit_reached=bound,
            projection_center_offsets_xy_px=fit.x[3:].tolist(),validation_split='spatial_96px_reference_cells',
            shared_geometry_training_exposure=n in calibration_exposures)
        quality['numerical_pass']=bool(shared_ok and fit.success and not bound and training_inliers.sum()>=8
            and confirmed.sum()>=4 and confirmed.mean()>=.6 and np.median(validation)<2
            and len(np.unique(sectors))>=3 and span>.5*radius)
        mapping=dict(kind='equisolid_shared_tangential_geometry',source_center_xy_px=center.tolist(),
            target_projection_center_xy_px=(center+fit.x[3:]).tolist(),radius_px=radius,
            ray_rotation_source_to_reference=Rotation.from_rotvec(fit.x[:3]).as_matrix().tolist(),
            shared_tangential_p1_p2=coefficients.tolist(),fitted_image_geometry_only=True)
        fitted[n]=(mapping,quality,errors)
    return fitted,shared

