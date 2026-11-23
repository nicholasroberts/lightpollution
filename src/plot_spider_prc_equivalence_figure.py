"""Publication figure drafts from saved PRC-equivalent spectral scenarios.

Summary uses saved numbers only. Optional matched-field maps require the
checked lock, fit arrays and measured Stokes on EcoVis; no fits are changed.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS={'HPS':'#D18A00','LED':'#0072B2'}
LABELS={'WClawn':'Worcester College lawn','oxfordbussinesspark':'Oxford Business Park'}


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def response_curve(report,lamp,energy,response='RH1_510'):
    fit=report['results'][response]['fits'][lamp]
    if fit['status']!='fitted' or fit['energy'] is None or fit['energy']<=0:
        raise ValueError('Figure requires a positive fitted light amount')
    k0=np.tanh(fit['natural_prc']/2);kt=np.tanh(report['target_maximum_absolute_prc']/2)
    dilution=(k0/kt-1)/fit['energy']
    return 2*np.arctanh(k0/(1+dilution*np.asarray(energy)))


def style(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(direction='out',length=3)


def letter(ax,value):ax.text(-.12,1.07,value,transform=ax.transAxes,fontweight='bold',fontsize=12)


def save(fig,output,name):
    for ext in ['png','pdf','svg']:
        fig.savefig(output/(name+'.'+ext),dpi=300,facecolor='white',bbox_inches='tight')
    plt.close(fig)


def plot_curve(ax,report,limit):
    energy=np.linspace(0,limit,700);target=report['target_maximum_absolute_prc']
    for lamp in ['HPS','LED']:
        amount=report['results']['RH1_510']['fits'][lamp]['energy']
        ax.plot(energy,response_curve(report,lamp,energy),color=COLORS[lamp],lw=2,label='Sodium (HPS)' if lamp=='HPS' else 'LED')
        ax.plot([amount,amount],[0,target],ls=':',lw=1,color=COLORS[lamp])
        ax.plot(amount,target,'o',color=COLORS[lamp],ms=5)
        ax.annotate(f'{amount:.2f}',(amount,target),xytext=(0,-18),textcoords='offset points',ha='center',color=COLORS[lamp],fontsize=9)
    ax.axhline(target,color='#333333',ls='--',lw=1)
    ax.text(.98,target+.016,f'Measured target {target:.3f}',ha='right',va='bottom',transform=ax.get_yaxis_transform(),fontsize=9)
    ax.set(xlim=(0,limit),ylim=(0,.7),ylabel='Maximum |PRC|')
    ax.set_xlabel('Added sky energy / lunar zenith energy\n352–700 nm',labelpad=8)
    ax.set_title(LABELS.get(report['site'],report['site']),loc='left',fontsize=11,pad=12)
    style(ax)


def summary(reports,a):
    fig,axes=plt.subplots(2,2,figsize=(10,7.6));fig.subplots_adjust(hspace=.75,wspace=.35,bottom=.15,top=.9)
    limit=max(r['results']['RH1_510']['fits']['HPS']['energy'] for r in reports)*1.2
    for ax,report,l in zip(axes[0],reports,['A','B']):plot_curve(ax,report,limit);letter(ax,l)
    axes[0,0].legend(frameon=False,fontsize=9,loc='upper right')
    first=reports[0];ax=axes[1,0]
    ratios=[first['results'][key]['led_to_hps_energy_ratio'] for key in ['RH1_510','D850_green_sensitivity']]
    for i,ratio in enumerate(ratios):
        ax.barh(i+.16,100,height=.28,color=COLORS['HPS'])
        ax.barh(i-.16,100*ratio,height=.28,color=COLORS['LED'])
        ax.text(100*ratio+2,i-.16,f'{100*ratio:.1f}',va='center',fontsize=9)
        ax.text(102,i+.16,'100',va='center',fontsize=9)
    ax.set_yticks([0,1]);ax.set_yticklabels(['RH1 model\n510 nm','D850 green\nsensitivity check'])
    ax.invert_yaxis();ax.set(xlim=(0,115),xlabel='Energy for the same PRC\n(sodium = 100)');style(ax);letter(ax,'C')
    ax=axes[1,1]
    from spider_spectral_lunar import pigment_proxy
    d=np.genfromtxt(a.camera_response,delimiter=',',names=True);w=d['wavelength_nm'];green=d['green_photon_response'];green/=green.max()
    ax.plot(w,pigment_proxy(w),color='#009E73',lw=2,label='RH1 510 nm model')
    ax.plot(w,green,color='#666666',lw=1.7,ls='--',label='Measured D850 green')
    ax.set(xlim=(350,700),ylim=(0,1.08),xlabel='Wavelength (nm)',ylabel='Relative photon sensitivity')
    ax.legend(frameon=False,fontsize=8,loc='upper right');style(ax);letter(ax,'D')
    reduction=first['results']['RH1_510']['led_energy_reduction_percent']
    axes[1,0].set_title(f'RH1: {reduction:.1f}% less LED energy',loc='left',fontsize=11,pad=12)
    fig.suptitle('Same receptor contrast under sodium and LED skyglow',fontsize=14,y=.985)
    fig.text(.5,.025,'Spectral scenarios: fixed Sun–Moon pattern + uniform unpolarized light.\nGreen-camera PRC used as the RH1 target approximation; the shared percentage follows from spectral weighting.',ha='center',fontsize=9,color='#444444')
    save(fig,a.output,'prc_equivalence_summary')


def block_mean(v,n):
    h,w=v.shape;v=v[:h//n*n,:w//n*n].reshape(h//n,n,w//n,n)
    count=np.isfinite(v).sum(axis=(1,3));total=np.nansum(v,axis=(1,3))
    return np.divide(total,count,out=np.full(total.shape,np.nan),where=count>0)


def maps(reports,a):
    locks=json.loads(a.reference_lock.read_text())['sites']
    checks=json.loads((a.reference_lock.parent/'reference_sha256.json').read_text())
    fig=plt.figure(figsize=(13.5,7));gs=fig.add_gridspec(2,5,width_ratios=[1,1,1,1.25,1.25],hspace=.85,wspace=.5)
    diagnostics=[];image=None;difference=None
    limit=max(r['results']['RH1_510']['fits']['HPS']['energy'] for r in reports)*1.2
    for row,report in enumerate(reports):
        site=report['site'];lock=locks[site];fitdir=a.fit_root/site/'fit';prcdir=a.prc_root/site
        fpath=fitdir/'twilight_moon_fit.json';apath=fitdir/'twilight_moon_fit_arrays.npz';tpath=prcdir/'prc_comparison.json';spath=prcdir/'prc_support.npz'
        for path,sha in [(fpath,report['source_fit_sha256']),(apath,report['source_fit_arrays_sha256']),
                         (tpath,report['source_prc_sha256']),(spath,report['source_support_sha256'])]:
            if digest(path)!=sha:raise ValueError('Source changed: '+str(path))
        fit=json.loads(fpath.read_text());target=json.loads(tpath.read_text())
        if digest(a.reference_lock)!=target['source_lock_sha256']:raise ValueError('Reference lock differs')
        stokespath=a.reference_lock.parent/site/'registration/linear_stokes.npz'
        if digest(stokespath)!=target['source_stokes_sha256'] or checks.get(str(stokespath.relative_to(a.reference_lock.parent)))!=digest(stokespath):raise ValueError('Measured Stokes checksum differs')
        with np.load(stokespath,allow_pickle=False) as z:I,Q,U=[block_mean(z[k],target['bin_size']) for k in ['I','Q','U']]
        with np.load(spath,allow_pickle=False) as z:mask=np.array(z['common_central_field'])
        with np.load(apath,allow_pickle=False) as z:
            natural_i=z['mixed_I']-fit['models']['mixed']['parameters'][1];pol=z['mixed_Q']+1j*z['mixed_U']
            if not np.allclose(np.hypot(Q,U)/I,z['measured_dolp'],equal_nan=True):raise ValueError('Measured arrays differ')
        result=report['results']['RH1_510'];cn=result['natural_response_per_band_energy']
        predicted={lamp:cn*abs(pol)/(cn*natural_i+result['fits'][lamp]['energy']*result['lamp_response_per_band_energy'][lamp]) for lamp in ['HPS','LED']}
        if not np.allclose(predicted['HPS'][mask],predicted['LED'][mask],rtol=1e-10,atol=1e-12):raise ValueError('Expected equivalent maps do not agree')
        observed=np.hypot(Q,U)/I;delta=predicted['HPS']-observed
        n=target['bin_size'];y,x=np.indices(I.shape);xy=np.stack([x*n+(n-1)/2,y*n+(n-1)/2],axis=-1)
        c=np.asarray(lock['reference_center_xy_px']);r=lock['eye_radius_px'];transform=np.asarray(fit['display_camera_to_north_xy'])
        north=xy@transform[:,:2].T+transform[:,2];xx=-(north[...,0]-c[0])/r;yy=-(north[...,1]-c[1])/r
        for col,(title,values) in enumerate([('Measured green\nRH1 approximation',observed),('PRC-matched RH1 model\nHPS and LED coincide',predicted['HPS']),('Model − measurement\nDoLP residual',delta)]):
            ax=fig.add_subplot(gs[row,col]);im=ax.pcolormesh(xx,yy,np.where(mask,values,np.nan),shading='nearest',rasterized=True,
                cmap='viridis' if col<2 else 'RdBu_r',vmin=0 if col<2 else -.08,vmax=.25 if col<2 else .08)
            ax.set(aspect='equal',xlim=(-1.1,1.1),ylim=(-1.1,1.1));ax.axis('off');ax.set_title(title,fontsize=9)
            ax.text(.02,.94,'N',transform=ax.transAxes,fontsize=8);ax.text(.01,.5,'E',transform=ax.transAxes,fontsize=8)
            letter(ax,chr(65+row*4+col))
            if col==0:ax.text(0,-1.38,LABELS.get(site,site)+'\nZenith-centred 60° field',ha='center',fontsize=9)
            if col<2:image=im
            else:difference=im
        ax=fig.add_subplot(gs[row,3:]);plot_curve(ax,report,limit);letter(ax,chr(68+row*4))
        if row==0:ax.legend(frameon=False,fontsize=8,loc='upper right')
        diagnostics.append(dict(site=site,dolp_rmse=float(np.sqrt(np.mean(delta[mask]**2))),
            dolp_mean_bias=float(np.mean(delta[mask])),maximum_hps_led_map_difference=float(np.max(abs(predicted['HPS'][mask]-predicted['LED'][mask])))))
    fig.subplots_adjust(bottom=.24,top=.9)
    cax=fig.add_axes([.12,.065,.21,.018]);fig.colorbar(image,cax=cax,orientation='horizontal',label='DoLP (common 0–0.25 scale)')
    cax=fig.add_axes([.39,.065,.12,.018]);fig.colorbar(difference,cax=cax,orientation='horizontal',label='DoLP difference (±0.08)')
    fig.text(.66,.035,'Sky view: north up, east left.\nOnly common measured field shown; no extrapolation.\nPRC matched; spatial DoLP is not forced to agree.',fontsize=8)
    fig.suptitle('Same receptor contrast, different required light amounts',fontsize=14)
    save(fig,a.output,'prc_equivalence_spatial')
    return diagnostics


def run(a):
    if len(a.sites)!=2:raise ValueError('This paired-site layout requires exactly two sites')
    reports=[];sources={}
    for site in a.sites:
        path=a.results/site/'prc_lamp_equivalence.json';r=json.loads(path.read_text())
        if r['site']!=site:raise ValueError('Site mismatch')
        if digest(a.camera_response)!=r['camera_response_sha256']:raise ValueError('Camera response changed')
        for response in ['RH1_510','D850_green_sensitivity']:
            for lamp in ['HPS','LED']:
                v=r['results'][response]['fits'][lamp]
                if v['status']!='fitted' or not np.isclose(v['achieved_prc'],r['target_maximum_absolute_prc'],atol=1e-10):raise ValueError('Target not matched')
        reports.append(r);sources[site]=digest(path)
    for key in ['RH1_510','D850_green_sensitivity']:
        if not np.isclose(reports[0]['results'][key]['led_to_hps_energy_ratio'],reports[1]['results'][key]['led_to_hps_energy_ratio']):raise ValueError('Summary cannot pool different spectral ratios')
    if reports[0]['spectral_band_nm']!=reports[1]['spectral_band_nm']:raise ValueError('Bands differ')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    a.output.mkdir(parents=True)
    summary(reports,a)
    supplied=[a.reference_lock,a.fit_root,a.prc_root]
    if any(supplied) and not all(supplied):raise ValueError('Spatial figure requires --reference-lock, --fit-root and --prc-root together')
    diagnostics=maps(reports,a) if all(supplied) else []
    with (a.output/'figure_values.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['site','response','lamp','required_band_energy','target_maximum_prc','led_to_hps_energy_ratio'])
        for r in reports:
            for response,values in r['results'].items():
                for lamp in ['HPS','LED']:writer.writerow([r['site'],response,lamp,values['fits'][lamp]['energy'],r['target_maximum_absolute_prc'],values['led_to_hps_energy_ratio']])
    (a.output/'figure_provenance.json').write_text(json.dumps(dict(source_reports_sha256=sources,script_sha256=digest(__file__),
        spatial_diagnostics=diagnostics,fit_changed=False,percentage_interpretation='Conditional spectral-model equivalence, not two independent estimates of a site-specific effect',
        map_interpretation='HPS/LED RH1 maps coincide at matched PRC under uniform unpolarized addition. Measured maps remain camera-green approximations.'),indent=2)+'\n')
    print('Saved figure drafts:',a.output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--sites',nargs='+',default=['WClawn','oxfordbussinesspark'])
    p.add_argument('--camera-response',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/d850_elf_response.csv')
    p.add_argument('--reference-lock',type=Path);p.add_argument('--fit-root',type=Path);p.add_argument('--prc-root',type=Path)
    a=p.parse_args()
    try:run(a)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))


if __name__=='__main__':main()
