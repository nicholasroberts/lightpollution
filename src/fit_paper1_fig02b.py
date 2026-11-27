"""Figures2B/2C: independent HPS/LED fits to site spatial camera-green DoLP."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_paper1_fig02a import spectral_inputs,integral
from paper1_figure_archive import archive_figure,digest


def fit_amount(intensity,polarized,observed,lamp_overlap):
    """Independently minimize equal-pixel spatial DoLP squared residuals."""
    i,p,y=[np.asarray(v,dtype=float).ravel() for v in (intensity,polarized,observed)]
    if i.size==0 or i.shape!=p.shape or i.shape!=y.shape or not np.isfinite(i).all() or not np.isfinite(p).all() or np.any(i<=0) or np.any(p<0):
        raise ValueError('Invalid natural Stokes support')
    if not np.isfinite(y).all() or np.any((y<0)|(y>1)) or not np.isfinite(lamp_overlap) or lamp_overlap<=0:
        raise ValueError('Invalid observed DoLP or spectral overlap')
    if not np.any(y>0):
        raise ValueError('All observed DoLP is zero: no finite lamp amount is identifiable')
    def residual(a):return p/(i+a[0]*lamp_overlap)-y
    def jacobian(a):return (-p*lamp_overlap/(i+a[0]*lamp_overlap)**2)[:,None]
    scale=float(np.mean(i)/lamp_overlap)
    # Each lamp has its own optimizer, physical amount, initialization and result.
    result=least_squares(residual,[scale],jac=jacobian,bounds=(0.,np.inf),
                         x_scale=scale,ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=1000)
    if not result.success:raise ValueError('Lamp optimization failed: '+result.message)
    amount=float(result.x[0])
    zero_cost=float(np.mean(residual([0.])**2))
    fitted_cost=float(np.mean(residual([amount])**2))
    if zero_cost<=fitted_cost:
        amount=0.
    prediction=p/(i+amount*lamp_overlap)
    return dict(amount=amount,achieved_mean_dolp=float(prediction.mean()),
                objective_mean_squared_DoLP_residual=float(np.mean((prediction-y)**2)),
                zero_ALAN_mean_squared_DoLP_residual=zero_cost,
                optimizer='scipy.optimize.least_squares; separate call per lamp',
                optimizer_success=bool(result.success),optimizer_message=result.message,
                optimizer_nfev=int(result.nfev),boundary_zero=bool(amount==0.))



def main(expected_site='Cottage',panel_letter='B'):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fit',type=Path,default=Path('data/processed/paper1_fig01c_joint_separation_01')/expected_site/'fit')
    p.add_argument('--north',type=Path,default=Path('data/processed/paper1_cottage_north_01/north') if expected_site=='Cottage' else Path('data/processed/paper1_museum_north_01'))
    p.add_argument('--reference-lock',type=Path,default=Path('data/processed/paper1_cottage_north_01/registration/locked_site_references.json') if expected_site=='Cottage' else Path('data/processed/spider_registration_green_all_sites_01/runs/museum/locked_site_references.json'))
    p.add_argument('--spectra-root',type=Path,default=Path('data/processed/paper1_spectral_inputs_02'))
    p.add_argument('--camera-response',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/d850_elf_response.csv')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--figure-copy',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.figure_copy.exists():p.error('Use fresh output and figure-copy directories')
    rp=a.fit/'solar_lunar_fit.json';ap=a.fit/'solar_lunar_fit_arrays.npz'
    fit=json.loads(rp.read_text());site=fit['site']
    if site!=expected_site:p.error('Expected '+expected_site+' for Figure2'+panel_letter)
    for path,key in [(a.reference_lock,'source_lock_sha256'),(a.north/'north_stokes.npz','source_north_arrays_sha256'),(a.north/'north_rotation.json','source_north_report_sha256')]:
        if digest(path)!=fit[key]:p.error('Source differs from final Figure1 fit: '+str(path))
    if not fit.get('neutral_separation_fitted') or not fit.get('shared_solar_lunar_separation'):p.error('Requires final joint Figure1 natural fit')
    lock=json.loads(a.reference_lock.read_text())['sites'][site]
    if lock.get('camera_channel')!='green':p.error('Camera-green reconstruction required')
    with np.load(ap,allow_pickle=False) as z:arr={k:np.array(z[k]) for k in z.files}
    f=fit['solar_fraction']
    for k in ('I','Q','U'):
        if not np.allclose(arr['combined_'+k],(1-f)*arr['lunar_'+k]+f*arr['solar_'+k],equal_nan=True):p.error('Saved natural mixture differs: '+k)
    wave,spectra,spectral_report=spectral_inputs(a.spectra_root)
    response=np.genfromtxt(a.camera_response,delimiter=',',names=True)
    if not {'wavelength_nm','green_photon_response'}<=set(response.dtype.names or []):p.error('Invalid camera response CSV')
    cw=response['wavelength_nm'];cg=response['green_photon_response']
    if not np.isfinite(cw).all() or np.any(np.diff(cw)<=0) or not np.isfinite(cg).all() or np.any(cg<0) or cg.max()<=0:p.error('Invalid camera response values')
    if wave[0]<cw[0] or wave[-1]>cw[-1]:p.error('Camera response does not cover the2A spectral band; no extrapolation')
    green=np.interp(wave,cw,cg)/cg.max();weight=wave*green
    overlaps={k:float(integral(v*weight,wave)) for k,v in spectra.items()}
    if overlaps['Moon']<=0:p.error('Zero natural camera-response overlap')
    x,y=arr['x_px'],arr['y_px'];c=np.asarray(lock['reference_center_xy_px']);r=lock['rim_radius_px']
    ni,nq,nu=[arr['combined_'+k] for k in ('I','Q','U')]
    mi,mq,mu=[arr['measured_'+k] for k in ('I','Q','U')]
    with np.errstate(divide='ignore',invalid='ignore'):
        measured=np.hypot(mq,mu)/mi;natural=np.hypot(nq,nu)/ni
    field=np.hypot(x-c[0],y-c[1])<lock['eye_radius_px']
    support=field & (mi>0) & (ni>0) & np.isfinite(measured) & (measured>=0) & (measured<=1) & np.isfinite(natural) & (natural>=0) & (natural<=1)
    if support.sum()<20:p.error('Insufficient common60degree field support')
    target=float(measured[support].mean())
    # Same natural spectral template for Sun and Moon preserves the green fit weights.
    camera_i=ni*overlaps['Moon'];camera_q=nq*overlaps['Moon'];camera_u=nu*overlaps['Moon']
    polarized=np.hypot(camera_q,camera_u)
    fits={};maps={};saved=dict(x_px=x,y_px=y,common_central_field=support,measured_DoLP=measured,
        natural_I=ni,natural_Q=nq,natural_U=nu,camera_natural_I=camera_i,camera_natural_Q=camera_q,camera_natural_U=camera_u,
        wavelength_nm=wave,moon_spectrum=spectra['Moon'],HPS_spectrum=spectra['HPS'],LED_spectrum=spectra['LED'],green_photon_response=green)
    for lamp in ('HPS','LED'):
        try:result=fit_amount(camera_i[support],polarized[support],measured[support],overlaps[lamp])
        except ValueError as exc:p.error(str(exc))
        amount=result['amount'];ii=camera_i+amount*overlaps[lamp]
        with np.errstate(divide='ignore',invalid='ignore'):d=polarized/ii
        residual=d[support]-measured[support]
        result.update(camera_added_I=float(amount*overlaps[lamp]),
            spatial_DoLP_rms_residual=float(np.sqrt(np.mean(residual**2))),
            mean_DoLP_residual=float(np.mean(residual)),fitted_mean_spatial_sd=float(d[support].std(ddof=1)))
        fits[lamp]=result;maps[lamp]=d
        saved[lamp+'_camera_I']=ii;saved[lamp+'_camera_Q']=camera_q;saved[lamp+'_camera_U']=camera_u;saved[lamp+'_DoLP']=d
    ah,al=fits['HPS']['amount'],fits['LED']['amount'];ratio=al/ah if ah>0 else None
    difference=None if ratio is None else 100*(1-ratio)
    with np.load(a.north/'north_stokes.npz',allow_pickle=False) as z:moon=np.array(z['moon_xy'])
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,4,figsize=(13,3.9),layout='constrained')
    titles=['Measured green','Natural: final Figure1 mix',f"HPS added\nBand radiance {ah:.3g}",f"LED added\nBand radiance {al:.3g}"]
    disc=np.hypot(x-c[0],y-c[1])<.985*r
    for ax,v,title in zip(axes,[measured,natural,maps['HPS'],maps['LED']],titles):
        valid=disc & np.isfinite(v) & (v>=0) & (v<=1)
        im=ax.pcolormesh(x,y,np.where(valid,v,np.nan),shading='nearest',cmap='viridis',vmin=0,vmax=1,rasterized=True)
        ax.add_patch(plt.Circle(c,lock['eye_radius_px'],fill=False,color='white',lw=.8));ax.plot(*moon,'k+',ms=4)
        ax.set(xlim=(c[0]-1.07*r,c[0]+1.07*r),ylim=(c[1]+1.07*r,c[1]-1.07*r),aspect='equal');ax.axis('off');ax.set_title(title,fontsize=10)
        for tx,ty,t in [(.5,.99,'N'),(.01,.5,'E'),(.5,.01,'S'),(.99,.5,'W')]:ax.text(tx,ty,t,transform=ax.transAxes,ha='center',va='center',fontsize=8)
    fig.colorbar(im,ax=axes.tolist(),shrink=.75,label='Camera-green DoLP')
    fig.suptitle(f'{site} — unpolarized HPS and LED additions to the final natural pattern\\n'
                 'Amounts relative to model natural sky band radiance; spatial camera-green DoLP fits',fontsize=11)
    fig.text(.01,.99,panel_letter,ha='left',va='top',fontweight='bold',fontsize=18)
    report=dict(subfigure='2'+panel_letter,site=site,target_mean_camera_green_dolp=target,
        target_spatial_sd=float(measured[support].std(ddof=1)),common_pixels=int(support.sum()),
        natural_mean_dolp=float(natural[support].mean()),natural_spatial_sd=float(natural[support].std(ddof=1)),
        solar_fraction=f,lunar_fraction=1-f,neutral_separation_deg=fit['neutral_separation_deg'],fits=fits,
        LED_to_HPS_band_radiance_ratio=ratio,LED_band_radiance_reduction_percent=difference,
        band_radiance_definition='Added sky radiant energy over shared band divided by total natural solar+lunar sky band energy per direction (uniform I=1); not lux, electrical power or lunar-only intensity',
        objective='Independently minimize equal-pixel mean squared spatial camera-green DoLP residuals for each lamp on the same common zenith-centred60deg field; binned Figure1 Stokes; no threshold0.05 mask',
        spectral_overlaps=overlaps,spectral_inputs=spectral_report,camera_response=dict(path=str(a.camera_response),sha256=digest(a.camera_response),convention='Published ELF green photon response; energy spectra weighted by wavelength*response; added polarizer not calibrated'),
        source_fit_report_path=str(rp),source_fit_arrays_path=str(ap),source_reference_lock_path=str(a.reference_lock),
        source_fit_report_sha256=digest(rp),source_fit_arrays_sha256=digest(ap),source_lock_sha256=digest(a.reference_lock),script_sha256=digest(__file__),
        colormaps=dict(DoLP='viridis'),dolp_colour_normalization='linear0–1',atmospheric_depolarization_multiplier=1.,extra_atmospheric_depolarization=False,
        limitations=['Uniform unpolarized addition leaves Q/U and AoP unchanged.',
        'Fixed scalar spectral overlaps make HPS/LED alternative fits mathematically equivalent in camera space; an identical amount ratio across sites can persist despite independent optimizations.',
        'The full spatial residual vector is minimized; the fitted mean is not constrained to equal the observed mean.',
        'Experimental statistics here are on the common binned model grid, not replacement original-camera Figure1 statistics.',
        'Common natural spectral template for solar and lunar components; no separate twilight-colour retrieval.',
        'No PRC or animal-response fit here;2E must use these same fitted physical lamp amounts.',
        'Fitted added radiance is conditional on the unity-depolarization baseline, spectra and uniform-sky approximation.'])
    a.output.mkdir(parents=True)
    arrays=a.output/'lamp_fit_arrays.npz';np.savez_compressed(arrays,**saved);report['source_output_arrays_sha256']=digest(arrays);report['source_output_arrays_path']=str(arrays)
    panel=a.output/'panel';panel.mkdir()
    for ext in ('png','pdf','svg'):fig.savefig(panel/('fig02'+panel_letter.lower()+'_'+site+'.'+ext),dpi=300)
    plt.close(fig)
    for dest in (a.output,panel):(dest/'lamp_DoLP_fit.json').write_text(json.dumps(report,indent=2)+'\n')
    with (panel/'lamp_amounts.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['lamp','added_band_radiance_relative_to_total_natural','mean_DoLP','spatial_DoLP_RMS_residual'])
        for lamp,result in fits.items():writer.writerow([lamp,result['amount'],result['achieved_mean_dolp'],result['spatial_DoLP_rms_residual']])
    archive_figure(panel,a.figure_copy)
    print(json.dumps(dict(target=target,fits=fits,LED_band_radiance_reduction_percent=difference),indent=2));print('Saved',a.figure_copy)


if __name__=='__main__':main()
