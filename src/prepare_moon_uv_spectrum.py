"""Audit and join measured near-full lunar source spectra, 352–750 nm.

Relative spectral shape only; QASUME and air-LUSI are TOA irradiance,
not measured ground-sky radiance. No extrapolation below validated data.
"""
from pathlib import Path
import argparse,hashlib,json,urllib.request
import numpy as np
from scipy.integrate import trapezoid

QURL='https://zenodo.org/records/21280043/files/DAV_QASUME_ToA_LunarIrradiance_Apr25.xlsx?download=1'
AURL='https://media.githubusercontent.com/media/usnistgov/air-lusi/main/data/air_lusi_spectra.nc'


def normalize(w,v):
    select=(w>=450)&(w<=480)
    if select.sum()<3 or np.any(~np.isfinite(v)) or np.any(v<=0):raise ValueError('Invalid lunar spectrum')
    return v/(trapezoid(v[select],w[select])/(w[select][-1]-w[select][0]))


def bin_means(w,v,edges):
    return np.array([trapezoid(np.interp(np.linspace(a,b,101),w,v),np.linspace(a,b,101))/(b-a) for a,b in zip(edges[:-1],edges[1:])])


def build(qpath,apath,out):
    import openpyxl,h5py
    sheet=openpyxl.load_workbook(qpath,data_only=True).worksheets[0]
    # Author's QASUME_Irradiance sheet: accepted absolute spectra only.
    # Last numeric wavelength row (550) has uncertainties but missing irradiance.
    qw=np.array([sheet.cell(r,1).value for r in range(23,221)],float)
    if not np.array_equal(qw,np.arange(352.,550.)):raise ValueError('Unexpected QASUME layout')
    scans=[];times=[];uncertainties=[]
    for col in [47,52,57,62,67]:
        v=np.array([sheet.cell(r,col).value for r in range(23,221)],float)
        scans.append(normalize(qw,v));times.append(sheet.cell(20,col).value.isoformat())
        uncertainties.append(np.array([sheet.cell(r,col+1).value for r in range(23,221)],float))
    q=np.median(scans,axis=0)
    with h5py.File(apath) as f:
        aw=f['Wavelength'][()];phases=f['signed_phase'][()]
        index=int(np.argmin(abs(abs(phases)-22.3)))
        av=f['Irradiance'][()][index];at=f['Time_String'][()][index].decode()
        au=f['Total_rel_err'][()][index]
    if not (aw[0]<=382 and aw[-1]>=750):raise ValueError('Unexpected air-LUSI support')
    av=normalize(aw,av)
    edges=np.arange(400.,541.,10.);ratio=bin_means(qw,q,edges)/bin_means(aw,av,edges)
    med=float(np.median(abs(ratio-1)));maximum=float(np.max(abs(ratio-1)))
    # Coarse-shape acceptance gate; not a statistical validation of UV accuracy.
    if med>.05 or maximum>.15:raise ValueError(f'Lunar spectral overlap failed: median {med:.3g}, max {maximum:.3g}')
    wave=np.arange(352.,751.);qa=np.interp(wave,qw,q);aa=np.interp(wave,aw,av)
    weight=np.clip((wave-400)/20,0,1);joined=(1-weight)*qa+weight*aa
    spread=np.std(scans,axis=0)/q
    out.mkdir(parents=True)
    np.savetxt(out/'moon_source_352_750.csv',np.c_[wave,joined],delimiter=',',header='wavelength_nm,relative_energy',comments='')
    np.savetxt(out/'overlap_10nm.csv',np.c_[edges[:-1],edges[1:],ratio],delimiter=',',header='bin_start_nm,bin_end_nm,qasume_over_air_lusi',comments='')
    np.savetxt(out/'qasume_selected_shapes.csv',np.c_[qw,q,spread],delimiter=',',header='wavelength_nm,median_relative_energy,relative_between_scan_std',comments='')
    np.savetxt(out/'air_lusi_selected_shape.csv',np.c_[aw,av,au],delimiter=',',header='wavelength_nm,relative_energy,relative_uncertainty',comments='')
    audit=dict(status='PASSED_COARSE_SPECTRAL_SHAPE_CHECK_NOT_SITE_VALIDATION',band_nm=[352,750],normalization_band_nm=[450,480],blend_nm=[400,420],overlap_bin_width_nm=10,overlap_median_absolute_fractional_difference=med,overlap_max_absolute_fractional_difference=maximum,acceptance_limits=dict(median=.05,max=.15),qasume=dict(url=QURL,doi='10.5281/zenodo.21280043',sha256=hashlib.sha256(qpath.read_bytes()).hexdigest(),sheet=sheet.title,scan_start_times_as_published=times,selection='Five complete accepted spectra Apr10–11, near-full phase; median normalized shape, not absolute calibration',valid_irradiance_nm=[352,549],uv_352_380_between_scan_relative_std_median=float(np.median(spread[qw<=380])),uv_reported_expanded_uncertainty_percent_median=float(np.median(np.array(uncertainties)[:,qw<=380]))),air_lusi=dict(url=AURL,doi='10.18434/mds2-3397',sha256=hashlib.sha256(apath.read_bytes()).hexdigest(),time_utc=at,phase_angle_deg=float(phases[index]),valid_nm=[float(aw[0]),float(aw[-1])],normalization='Use standardized-distance Irradiance; no distance correction needed for relative shape; known metadata issue5 recorded'),limitations=['Different dates and near-full phases; not exact full moon and not site-specific spectral retrieval','TOA source irradiance; atmospheric extinction and scattering still required for ground sky','Overlap check assesses coarse visible shape only, not UV absolute accuracy','UV coverage ends at352nm; no extrapolation below352','QASUME uncertainty includes atmospheric-correction assumptions; between-scan spread is not a confidence interval'])
    (out/'moon_source_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    axes[0].plot(qw,q,label='QASUME Apr10–11 median');axes[0].plot(aw,av,label='air-LUSI '+at[:10]);axes[0].plot(wave,joined,'k--',label='Composite');axes[0].set(xlim=(350,750),xlabel='Wavelength (nm)',ylabel='Relative energy; mean450–480nm=1');axes[0].legend(fontsize=8)
    axes[1].plot((edges[:-1]+edges[1:])/2,100*(ratio-1),'o-');axes[1].axhline(0,color='gray');axes[1].set(xlabel='10nm bin centre (nm)',ylabel='QASUME / air-LUSI difference (%)',title='Coarse shape check; median abs '+f'{100*med:.2f}%')
    fig.savefig(out/'moon_spectrum_overlap.png',dpi=160);plt.close(fig)
    print(f'Lunar overlap: median {100*med:.2f}%; max {100*maximum:.2f}%; UV scan spread {100*audit["qasume"]["uv_352_380_between_scan_relative_std_median"]:.2f}%')
    return out/'moon_source_352_750.csv'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--qasume',type=Path);p.add_argument('--air-lusi',type=Path);a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        paths=[]
        for path,url,name in [(a.qasume,QURL,'q.xlsx'),(a.air_lusi,AURL,'a.nc')]:
            if path is None:path=Path(td)/name;path.write_bytes(urllib.request.urlopen(url,timeout=90).read())
            paths.append(path)
        build(*paths,a.output)

if __name__=='__main__':main()
