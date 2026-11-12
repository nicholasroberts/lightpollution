"""Independent radiance-proxy site contrast and DNB-overlap spectral scenarios.
No measured DoP, AoP or PRC read or fitted. Absolute ALAN remains a sensitivity range.
"""
from pathlib import Path
import argparse,csv,json,urllib.request,hashlib,io
import numpy as np
from scipy.integrate import trapezoid
URL='https://ncc.nesdis.noaa.gov/documents/documentation/GT_F1_RSR_DNBLGS_v5_response_orbit6557_ccm.dat'
ROOT=Path(__file__).resolve().parents[1]


def overlap_conversion(source,rsr,exponent=2.):
    """Sky450–750nm energy per response-weighted source energy within common support."""
    if np.any(source[:,1]<0) or np.any(np.diff(source[:,0])<=0):raise ValueError('Invalid source array')
    lo=max(source[0,0],rsr[0,0]);hi=min(source[-1,0],rsr[-1,0])
    wave=np.unique(np.r_[np.arange(lo,hi,1.),hi,source[(source[:,0]>=lo)&(source[:,0]<=hi),0],rsr[(rsr[:,0]>=lo)&(rsr[:,0]<=hi),0]])
    signal=np.interp(wave,source[:,0],source[:,1]);response=np.maximum(np.interp(wave,rsr[:,0],rsr[:,1]),0)
    denom=trapezoid(signal*response,wave)
    visible=np.arange(450.,751.);numer=trapezoid(np.interp(visible,source[:,0],source[:,1])*(550/visible)**exponent,visible)
    if denom<=0:raise ValueError('No DNB overlap')
    return float(numer/denom),[float(lo),float(hi)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectra',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--references',type=Path,default=ROOT/'data/reference/prc_spider_18_sites_verified.csv')
    p.add_argument('--rsr-file',type=Path);p.add_argument('--skyglow-ratio',type=float,help='Independent ground-sky ratio override; otherwise use radiance ratio as conditional proxy')
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    rows={r['site_id']:r for r in csv.DictReader(a.references.open())};c=float(rows['s14']['original_multiyear_monthly_radiance_mean']);o=float(rows['s02']['original_multiyear_monthly_radiance_mean'])
    ratio=o/c if a.skyglow_ratio is None else a.skyglow_ratio
    if not np.isfinite(ratio) or ratio<=0:p.error('Positive site contrast required')
    raw=a.rsr_file.read_bytes() if a.rsr_file else urllib.request.urlopen(URL,timeout=60).read()
    rsr=np.loadtxt(io.BytesIO(raw));sp=json.loads((a.spectra/'spectral_provenance.json').read_text())
    expo=sp['alan_transfer_exponent'];conversions={};support={}
    for label in ['HPS','LED_mid']:
        source=np.genfromtxt(a.spectra/(label+'_author_source.csv'),delimiter=',',names=True)
        conversions[label],support[label]=overlap_conversion(np.column_stack([source[n] for n in source.dtype.names]),rsr,expo)
    a.output.mkdir(parents=True)
    np.savetxt(a.output/'noaa_snpp_dnb_rsr.csv',rsr,delimiter=',',header='wavelength_nm,relative_response_as_published',comments='')
    for label in conversions:
        payload=dict(status='INDEPENDENT_FORWARD_SCENARIOS_NOT_DOP_FIT',normalization_site='Cottage',normalization_meaning='Scenario level = Oxford HPS sky450–750nm ALAN divided by Cottage capture zenith natural intensity; common reference for both sites/spectra',alan_scale_by_site={'Cottage':1/ratio,'oxfordbussinesspark':1.},spectral_band_factor=conversions[label]/conversions['HPS'],source_spectrum=label,site_contrast=dict(oxford_over_cottage=ratio,basis='multiyear upward radiance proxy, not measured skyglow' if a.skyglow_ratio is None else 'user-supplied independent ground-sky ratio',cottage_radiance=c,oxford_radiance=o),dnb=dict(source_url=URL,sha256=hashlib.sha256(raw).hexdigest(),overlap_nm=support[label],sky_energy_per_dnb_weighted_source=conversions[label],negative_response_samples_zeroed_for_integration=int((rsr[:,1]<0).sum())),recommended_scenario_levels=[.3,1.,3.],limitations=['Absolute upward-to-ground sky transfer uncalibrated; ALAN levels are sensitivity scenarios, not measured site intensity','Radiance ratio does not include spatial convolution, shielding, aerosol differences or different historical observing times','DNB-overlap normalization truncates at source support900nm; not full operational VIIRS calibration','Same assumed wavelength transfer and atmospheric geometry for HPS/LED; no DoP-targeted normalization','Same spectral hypothesis applied at both sites; real lamp mixtures unknown'])
        (a.output/(label+'_site_alan.json')).write_text(json.dumps(payload,indent=2)+'\n')
    print(f'Oxford/Cottage radiance-proxy contrast {ratio:.3f}; LED/HPS sky-band factor at equal DNB overlap {conversions["LED_mid"]/conversions["HPS"]:.3f}')

if __name__=='__main__':main()
