"""Audit actual Figure 2 spectra and camera weighting; no fitting or artwork changes."""
import argparse
import json
from pathlib import Path
import numpy as np
from plot_paper1_fig02a import spectral_inputs,integral
from paper1_figure_archive import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectra-root',type=Path,required=True)
    p.add_argument('--camera-response',type=Path,default=Path(__file__).resolve().parents[1]/'data/reference/d850_elf_response.csv')
    p.add_argument('--fit-results',type=Path,nargs='*',default=[])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output file')
    wave,spectra,provenance=spectral_inputs(a.spectra_root)
    camera=np.genfromtxt(a.camera_response,delimiter=',',names=True)
    cw,cg=camera['wavelength_nm'],camera['green_photon_response']
    if not np.isfinite(cw).all() or np.any(np.diff(cw)<=0) or not np.isfinite(cg).all() or np.any(cg<0) or cg.max()<=0:
        p.error('Invalid camera response')
    if wave[0]<cw[0] or wave[-1]>cw[-1]:p.error('Spectra exceed camera response support')
    green=np.interp(wave,cw,cg)/cg.max()
    overlaps={k:float(integral(v*wave*green,wave)) for k,v in spectra.items()}
    if any(v<=0 for v in overlaps.values()):p.error('Nonpositive spectral overlap')
    checks=[]
    for root in a.fit_results:
        path=root/'lamp_DoLP_fit.json';fit=json.loads(path.read_text())
        same_input=fit['spectral_inputs']['source_provenance_sha256']==provenance['source_provenance_sha256']
        same_files=all(fit['spectral_inputs']['source_files'][k]['sha256']==provenance['source_files'][k]['sha256'] for k in ('HPS','LED'))
        same_camera=fit['camera_response']['sha256']==digest(a.camera_response)
        same_overlap=all(np.isclose(fit['spectral_overlaps'][k],v,rtol=1e-10,atol=0) for k,v in overlaps.items())
        ok=bool(same_input and same_files and same_camera and same_overlap)
        checks.append(dict(site=fit['site'],fit_path=str(path),input_hashes_and_overlaps_match=ok))
        if not ok:p.error('Saved fit inputs/overlaps differ from selected audit inputs: '+str(path))
    report=dict(status='NUMERICAL_INPUT_CHECK_PASSED_NOT_SITE_SPECTRAL_VALIDATION',
        spectral_inputs=provenance,camera_response_path=str(a.camera_response),
        camera_response_sha256=digest(a.camera_response),
        camera_convention='Counts per incident photon: multiply spectral energy by wavelength/(hc) once; hc cancels in relative Stokes ratios.',
        camera_overlap_relative_scale=overlaps,fit_checks=checks,
        limitations=['Generic emission templates plus declared approximate sky transfer, not measured pure site skyglow.',
                     'Read inherited provenance for actual generated band and transfer exponent; current script defaults do not establish historical file settings.',
                     'Added polarizer spectral transmission and modulation efficiency are not included in the ELF calibration.',
                     'Identical natural spectra for lunar and solar components remain an approximation.',
                     'No absolute emitted flux, lamp mixture, electrical power or site spectral validation is recovered.'])
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(band_nm=provenance['spectral_band_nm'],
        inherited_settings=provenance['inherited_spectral_provenance'],
        fit_checks=checks,output=str(a.output)),indent=2))

if __name__=='__main__':main()
