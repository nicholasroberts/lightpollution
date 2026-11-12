"""Retrieve author numerical spectra; prepare a VISIBLE-ONLY spectral pilot.

No UV extrapolation. NIST2013 Table1 moon irradiance is a source spectrum,
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
    p.add_argument('--alan-transfer-exponent',type=float,default=2.)
    a=p.parse_args()
    if a.output.exists():p.error('Use fresh output directory')
    if not np.isfinite(a.alan_transfer_exponent) or not 0<=a.alan_transfer_exponent<=4:p.error('Transfer exponent must be0..4')
    raw=a.archive.read_bytes() if a.archive else urllib.request.urlopen(URL,timeout=60).read()
    z=zipfile.ZipFile(io.BytesIO(raw));out=a.output;out.mkdir(parents=True)
    wave=np.arange(450.,751.);natural=np.interp(wave,NIST[:,0],NIST[:,1])*(550/wave)**4
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
    uvwave=np.arange(300.,751.)
    for center in [365.,380.]:
        for amplitude in [0.,.29,1.]:
            curve=pigment_proxy(uvwave,uv_peak_nm=center,uv_amplitude=amplitude)
            np.savetxt(out/f'proxy_510_uv{center:g}_amplitude{amplitude:g}.csv',np.column_stack([uvwave,curve]),delimiter=',',header='wavelength_nm,relative_photon_sensitivity',comments='')
    payload=dict(status='VISIBLE_ONLY_SOURCE_TRANSFER_SCENARIOS_NOT_MEASURED_PURE_SKYGLOW',source_url=URL,archive_sha256=hashlib.sha256(raw).hexdigest(),source_paper_doi='10.1093/mnras/stz3260',moon_reference='NIST2013 doi10.6028/jres.118.020 Table1',sources=sources,pilot_band_nm=[450,750],natural_transfer_exponent=4,alan_transfer_exponent=a.alan_transfer_exponent,pigment=dict(alpha_peak_nm=510,uv_peak_nm=365,uv_amplitude=.29,uv_width_nm=40,uv_amplitude_status='illustrative_not_measured',template_doi='10.1007/s00359-010-0568-7'),limitations=['No natural UV spectrum ingested; spectral PRC predictions restricted450–750nm','Generic lamp spectra are source emission, not measured pure skyglow','Power-law transfers are sensitivity scenarios, not radiative-transfer retrieval','NIST table interpolated coarsely; source near full moon2012-11-30, not exact full phase or site atmosphere','All inputs relative energy except pigment photon sensitivity; lambda factor converts energy weighting to photons','No inference of site lamp mixture, CCT or ALAN amount'])
    (out/'spectral_provenance.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(f'Wrote visible pilot inputs, measured sky examples and UV proxy curves to {out}')

if __name__=='__main__':main()
