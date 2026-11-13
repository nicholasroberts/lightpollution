"""Prepare measured lunar352–750nm source and lamp spectral scenarios.

No UV extrapolation. Composite QASUME/air-LUSI moon irradiance is a source spectrum,
not diffuse sky radiance. Rayleigh lambda^-4 natural transfer and optional
ALAN lambda^-exponent are declared scenario approximations, not measurements.
"""
from pathlib import Path
import argparse,io,json,zipfile,urllib.request,hashlib
import numpy as np
from spider_spectral_lunar import pigment_proxy

URL='https://davinci.fmph.uniba.sk/~kundracik1/SpectrumAnalyser/SpectrumAnalyzer_v1.1_2025_02_27_incl_Windows_executables.zip'
NIST=np.array([[449.7,2.348,.85],[499.9,2.395,.56],[550.,2.633,.45],[600.2,2.669,.44],[650.1,2.598,.40],[702.8,2.474,.38],[750.,2.314,.37],[850.2,1.870,.36],[1000.2,1.387,.54]])
LAMPS={'HPS':'HPS_High_Pressure_Sodium_GENERIC.txt','LED_warm':'LED_1833K-2877K_GENERIC.txt','LED_mid':'LED_2958K-3969K_GENERIC.txt','LED_cool':'LED_4225K-9089K_GENERIC.txt'}


def numeric_rows(raw):
    rows=[]
    for line in raw.decode('utf-8',errors='replace').splitlines():
        if line.lstrip().startswith('#'):continue
        try:row=[float(x) for x in line.split()[:2]]
        except ValueError:continue
        if len(row)==2:rows.append(row)
    a=np.asarray(rows)
    if a.ndim!=2 or len(a)<2 or np.any(~np.isfinite(a)) or np.any(np.diff(a[:,0])<=0):raise ValueError('Invalid author spectrum')
    return a


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--archive',type=Path,help='Previously downloaded author archive, otherwise retrieve HTTPS')
    p.add_argument('--moon-source',type=Path,help='Audited source CSV; otherwise download and audit QASUME/air-LUSI')
    p.add_argument('--alan-transfer-exponent',type=float,default=2.)
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    if not np.isfinite(a.alan_transfer_exponent) or not 0<=a.alan_transfer_exponent<=4:p.error('Transfer exponent must be0..4')
    raw=a.archive.read_bytes() if a.archive else urllib.request.urlopen(URL,timeout=60).read()
    z=zipfile.ZipFile(io.BytesIO(raw));out=a.output;out.mkdir(parents=True)
    if a.moon_source:
        moon_path=a.moon_source
    else:
        import subprocess,sys
        subprocess.run([sys.executable,str(Path(__file__).with_name('prepare_moon_uv_spectrum.py')),'--output',str(out/'moon_audit')],check=True)
        moon_path=out/'moon_audit/moon_source_352_750.csv'
    moon=np.genfromtxt(moon_path,delimiter=',',names=True)
    wave=np.arange(352.,751.)
    if moon['wavelength_nm'][0]>wave[0] or moon['wavelength_nm'][-1]<wave[-1]:raise ValueError('Moon source must cover352–750nm without extrapolation')
    natural=np.interp(wave,moon['wavelength_nm'],moon['relative_energy'])*(550/wave)**4
    np.savetxt(out/'moon_source_352_750.csv',np.column_stack([moon['wavelength_nm'],moon['relative_energy']]),delimiter=',',header='wavelength_nm,relative_energy',comments='')
    audit_path=moon_path.with_name('moon_source_audit.json')
    if not audit_path.exists():raise ValueError('Moon source audit JSON required beside CSV')
    moon_audit=json.loads(audit_path.read_text())
    if moon_audit['status']!='PASSED_COARSE_SPECTRAL_SHAPE_CHECK_NOT_SITE_VALIDATION':raise ValueError('Moon audit not passed')
    np.savetxt(out/'nist2013_table1.csv',NIST,delimiter=',',header='wavelength_nm,toa_irradiance_uW_m2_nm,uncertainty_percent',comments='')
    # Photon sensitivity: multiply energy radiance by lambda; common hc cancels.
    sensitivity=wave*pigment_proxy(wave)
    sources={}
    for label,name in LAMPS.items():
        member='SpectrumAnalyzer_v1.1/kernel/spectra_original/'+name
        source=numeric_rows(z.read(member))
        if np.any(source[:,1]<0):raise ValueError('Negative lamp data: inspect before using')
        np.savetxt(out/(label+'_author_source.csv'),source,delimiter=',',header='wavelength_nm,relative_energy',comments='')
        alan=np.interp(wave,source[:,0],source[:,1])*(550/wave)**a.alan_transfer_exponent
        np.savetxt(out/(label+'_visible_input.csv'),np.column_stack([wave,natural,alan,sensitivity]),delimiter=',',header='wavelength_nm,natural,alan,sensitivity',comments='')
        sources[label]=dict(archive_member=member,sha256=hashlib.sha256(z.read(member)).hexdigest(),range_nm=source[[0,-1],0].tolist())
    for name in ['bratislava_2010.txt','vienna_2012.txt','vienna_2018.txt']:
        member='SpectrumAnalyzer_v1.1/kernel/experimental_data/'+name;source=numeric_rows(z.read(member))
        np.savetxt(out/(name.removesuffix('.txt')+'_measured_sky.csv'),source,delimiter=',',header='wavelength_nm,relative_measured_signal',comments='')
        sources[name]=dict(archive_member=member,negative_samples=int((source[:,1]<0).sum()),range_nm=source[[0,-1],0].tolist(),note='Mixed-source observed sky; negative noise retained; not used as pure HPS/LED')
    np.savetxt(out/'pigment_510_alpha.csv',np.column_stack([wave,pigment_proxy(wave)]),delimiter=',',header='wavelength_nm,relative_photon_sensitivity',comments='')
    payload=dict(status='UV_VISIBLE_SOURCE_TRANSFER_SCENARIOS_NOT_MEASURED_PURE_SKYGLOW',source_url=URL,archive_sha256=hashlib.sha256(raw).hexdigest(),source_paper_doi='10.1093/mnras/stz3260',moon_reference='QASUME + air-LUSI measured near-full source; audited composite',moon_audit=moon_audit,sources=sources,pilot_band_nm=[352,750],natural_transfer_exponent=4,alan_transfer_exponent=a.alan_transfer_exponent,pigment=dict(alpha_peak_nm=510,uv_amplitude=0,template='Stavenga green alpha band only; phenomenological RH1-type proxy, molecular identity unassigned',template_doi='10.1007/s00359-010-0568-7'),limitations=['No extrapolation below352nm; no added UV receptor or Gaussian peak','Generic lamp spectra are source emission, not measured pure skyglow; weak UV lamp samples may include noise','Power-law transfers are sensitivity scenarios, not radiative-transfer retrieval; UV atmospheric attenuation not independently calibrated','Composite source uses different near-full dates; not an exact full-moon or site spectral measurement','All inputs relative energy except pigment photon sensitivity; lambda factor converts energy weighting to photons','Legacy filenames ending visible_input.csv retained for workflow compatibility; actual band is in metadata','No inference of site lamp mixture, CCT or ALAN amount'])
    (out/'spectral_provenance.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(f'Wrote352–750nm inputs with single510nm alpha pigment to {out}')

if __name__=='__main__':main()
