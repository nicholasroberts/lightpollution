#!/usr/bin/env python3
"""Independent moth DRA-PRC SSasymp calibration parallel to the locked spider fit.

17 original R vector observations align with R list.files() minus positions
1 (Bristol.csv) and 5 (Gavarnie.csv); both are excluded BEFORE alignment.
Uses the 18-site, independently geolocated and dated annual EOG values recovered
for the spider analysis, but fits the MOTH response values separately.

Slabhouse (raw annual EOG 0 both 2018 and 2019, original monthly >0) remains
in archival data and full-17 sensitivity, but is excluded from the preferred
observation-year moth model for the same source-QA reason as in the spider fit.
No biological transfer coefficients are shared with the spider model.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FILENAME_ORDER = [
    'ChurchHill.csv', 'Cottage.csv', 'EastHorrington.csv',
    'minifactory.csv', 'museum.csv', 'opensquare.csv', 'Ounsdale.csv',
    'oxfordbussinesspark.csv', 'PennCC.csv', 'Skillinge.csv',
    'slabhouse.csv', 'southhorrington.csv', 'stannes.csv', 'Waitrose.csv',
    'WClawn.csv', 'WCpitch.csv', 'zoology.csv',
]
ORIGINAL_PRCS = np.array([
    .407, .603, .692, .239, .159, .093, .448, .206, .378, .539,
    .620, .629, .142, .211, .071, .109, .216,
])
SLABHOUSE = 'slabhouse.csv'


def project_path(s):
    p = Path(s)
    return p if p.is_absolute() else ROOT/p


def predict(x, z):
    # Exactly R SSasymp: yf + (y0-yf)*exp(-exp(log_alpha)*x).
    yf, y0, log_alpha = z
    return yf+(y0-yf)*np.exp(-np.exp(log_alpha)*np.asarray(x,float))


def fit(x, y):
    x, y = np.asarray(x,float), np.asarray(y,float)
    if len(x)<4 or not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError('SSasymp requires at least 4 complete observations')
    # Match the 9 multi-start fits used by src/fit_prc_observation_year.py.
    attempts = [least_squares(lambda z: predict(x,z)-y,
                 [yf, max(y), np.log(k)],max_nfev=4000)
                for yf in (.10,.15,.20) for k in (.04,.10,.25)]
    z = min(attempts, key=lambda a: float(np.dot(a.fun,a.fun))).x
    residual = y-predict(x,z)
    if not (0 < z[0] < z[1] < 1 and np.exp(z[2])>0):
        raise ValueError('Moth fit is not a bounded decreasing PRC response')
    return {'z':z, 'yf':float(z[0]),'y0':float(z[1]),
            'k':float(np.exp(z[2])),'R2':float(1-np.sum(residual**2)/np.sum((y-y.mean())**2)),
            'RMSE':float(np.sqrt(np.mean(residual**2))), 'n':len(x)}


def validate_inputs(moth, sites):
    req={'source_order','source_filename','moth_DRAc_original_R'}
    if not req.issubset(moth):
        raise ValueError('Moth R data columns missing')
    if len(moth)!=17 or moth.source_filename.tolist()!=FILENAME_ORDER or \
            not np.array_equal(moth.source_order.to_numpy(int),np.arange(1,18)):
        raise ValueError('Moth observations do not match original R ordering')
    if not np.allclose(moth.moth_DRAc_original_R, ORIGINAL_PRCS,rtol=0,atol=1e-12):
        raise ValueError('Moth DRAc values differ from original 17-entry R vector')
    required={'source_filename','site_id','field_date','field_year',
              'eog_corrected_field_year','original_multiyear_monthly_radiance_mean',
              'eog_2018_row','eog_2018_col','eog_2019_row','eog_2019_col'}
    if not required.issubset(sites):
        raise ValueError('Site-level EOG file is missing required fields')
    if len(sites)!=18 or not sites.source_filename.is_unique:
        raise ValueError('Expected 18 unique spider-derived source filename records')
    if sites.source_filename.isin(FILENAME_ORDER).sum()!=17:
        raise ValueError('One or more moth filenames missing from annual EOG sites')
    d=moth.merge(sites,on='source_filename',validate='one_to_one',how='left',sort=False)
    if d.site_id.isna().any() or len(d)!=17:
        raise ValueError('Unresolved moth site alignment')
    date_year=pd.to_datetime(d.field_date).dt.year
    if not date_year.eq(d.field_year).all() or not d.field_year.isin([2018,2019]).all():
        raise ValueError('Inconsistent or unknown field dates for moth site')
    if (int(d.field_year.eq(2018).sum()),int(d.field_year.eq(2019).sum()))!=(2,15):
        raise ValueError('Expected 2 moth observation sites in 2018 and 15 in 2019')
    if not np.isfinite(d.eog_corrected_field_year.to_numpy(float)).all():
        raise ValueError('Unusable field-year EOG radiance in moth input')
    d['EOG_pixel_group']=[f'{year}:{int(row)}:{int(col)}' for year,row,col in zip(
        d.field_year,
        np.where(d.field_year.eq(2018),d.eog_2018_row,d.eog_2019_row),
        np.where(d.field_year.eq(2018),d.eog_2018_col,d.eog_2019_col))]
    zero=d.loc[d.source_filename.eq(SLABHOUSE),'eog_corrected_field_year']
    if len(zero)!=1 or zero.iloc[0]!=0.0:
        raise ValueError('Slabhouse pixel QA discrepancy changed; reconsider exclusion')
    d['included_in_preferred_moth_fit']=d.source_filename.ne(SLABHOUSE)
    return d


def heldout_by_pixel(frame,xcol):
    errors=[]
    for key in frame.EOG_pixel_group.unique():
        test=frame.EOG_pixel_group.eq(key)
        trained=fit(frame.loc[~test,xcol],frame.loc[~test,'moth_DRAc_original_R'])
        estimated=predict(frame.loc[test,xcol],trained['z'])
        errors.extend((frame.loc[test,'moth_DRAc_original_R'].to_numpy()-estimated)**2)
    return float(np.sqrt(np.mean(errors)))


def self_test():
    assert len(FILENAME_ORDER)==17 and len(set(FILENAME_ORDER))==17
    assert FILENAME_ORDER[0]=='ChurchHill.csv' and FILENAME_ORDER[-1]=='zoology.csv'
    assert 'Gavarnie.csv' not in FILENAME_ORDER and 'Bristol.csv' not in FILENAME_ORDER
    x=np.array([0,1,5,10,20,40,60],float)
    true=np.array([.1,.7,np.log(.08)])
    info=fit(x,predict(x,true))
    assert np.allclose(info['z'],true,rtol=0,atol=1e-5)
    print('Moth R ordering and SSasymp synthetic data: PASS',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--moth-r',default='data/reference/prc_moth_17_observations_from_original_r.csv')
    ap.add_argument('--site-eog',default='data/processed/prc_spider_18_site_refit_confirmed/spider18_site_eog_corrected_2018_2019.csv')
    ap.add_argument('--spider-config',default='config/prc_transfer.yaml')
    ap.add_argument('--out',default='data/processed/prc_moth_observation_year')
    args=ap.parse_args()
    self_test()
    if args.self_test:return
    d=validate_inputs(pd.read_csv(project_path(args.moth_r)),
                      pd.read_csv(project_path(args.site_eog)))
    original='original_multiyear_monthly_radiance_mean'
    eog='eog_corrected_field_year'
    chosen=d[d.included_in_preferred_moth_fit].copy()
    assert len(chosen)==16 and chosen.EOG_pixel_group.nunique()==15
    specs=[('Original R monthly mean, all 17',d,original),
           ('Corrected EOG field year, all 17',d,eog),
           ('Original R monthly mean, excluding Slabhouse',chosen,original),
           ('Corrected EOG field year, excluding Slabhouse',chosen,eog)]
    fits={}; records=[]
    for name,frame,xx in specs:
        info=fit(frame[xx],frame.moth_DRAc_original_R)
        info['pixel_grouped_CV_RMSE']=heldout_by_pixel(frame,xx)
        records.append({'model':name,'x_column':xx,'distinct_satellite_pixels':
                        int(frame.EOG_pixel_group.nunique()),**{
                        k:v for k,v in info.items() if k!='z'}})
        fits[name]=info
    selected=fits[specs[-1][0]]
    d['moth_predicted_preferred_16']=predict(d[eog],selected['z'])
    d['moth_residual_preferred_16']=d.moth_DRAc_original_R-d.moth_predicted_preferred_16
    out=project_path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    d.to_csv(out/'moth_prc_site_alignment_and_residuals.csv',index=False)
    pd.DataFrame(records).to_csv(out/'moth_prc_fit_comparison.csv',index=False)

    import yaml
    with project_path(args.spider_config).open(encoding='utf-8') as f:
        spider=yaml.safe_load(f)
    if spider.get('status')!='locked' or spider['fit']['included_prc_observations']!=17:
        raise ValueError('Spider comparison expected locked 17-site transfer')
    sp=spider['parameters']
    fig,ax=plt.subplots(figsize=(10,6))
    xx=np.linspace(0,65,400)
    ax.scatter(d[eog],d.moth_DRAc_original_R,alpha=.55,marker='o',label='Moth: all original 17 observations')
    ax.scatter(chosen[eog],chosen.moth_DRAc_original_R,marker='o',facecolors='none',
               edgecolors='black',s=70,label='Moth: 16 observations retained')
    for name,kind in [(specs[0][0],'--'),(specs[-1][0],'-')]:
        info=fits[name]
        ax.plot(xx,predict(xx,info['z']),linestyle=kind,lw=2,
                label=f'Moth {name} (R²={info["R2"]:.3f})')
    spider_y=sp['yf']+(sp['y0']-sp['yf'])*np.exp(-sp['k']*xx)
    ax.plot(xx,spider_y,lw=2,linestyle=':',label='Spider: independently locked 17-site fit')
    ax.set(xlabel='VIIRS annual observation-year calibrated radiance (nW cm⁻² sr⁻¹)',
           ylabel='Polarization receptor contrast (PRC)',
           title='Original moth DRA PRC: independent SSasymp transfer')
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(out/'moth_prc_and_spider_fit_comparison.pdf',bbox_inches='tight')
    fig.savefig(out/'moth_prc_and_spider_fit_comparison.png',bbox_inches='tight',dpi=170)
    plt.close(fig)

    lines=[
      'INDEPENDENT MOTH PRC OBSERVATION-YEAR EOG CALIBRATION',
      'R original: list.files(path), filenames[-c(1,5)] = exclude Bristol and Gavarnie',
      'Original 17 DRAc PRCs preserved; 2 sites sampled in 2018, 15 in 2019.',
      'Same confirmed site coordinates, measurement years, EOG 2017+ calibration as spider.',
      'Slabhouse raw annual EOG pixel zero in both 2018 and 2019; exclude from',
      'preferred moth regression ONLY, preserve measured PRC=0.620 in reference.',
      'Preferred candidate n=16; 15 distinct EOG satellite pixels (Open Square/St Annes share).',
      'PRC(R)=yf+(y0-yf)*exp(-kR), fit exactly as spider R SSasymp.',
      f'Preferred yf={selected["yf"]:.12f}; y0={selected["y0"]:.12f}; k={selected["k"]:.12f}',
      f'Preferred PRC(R) = {selected["yf"]:.9f} + '
      f'{selected["y0"]-selected["yf"]:.9f} * exp(-{selected["k"]:.9f} * R)',
      f'Preferred R2={selected["R2"]:.7f}; RMSE={selected["RMSE"]:.7f}; '
      f'pixel-blocked CV RMSE={selected["pixel_grouped_CV_RMSE"]:.7f}',
      '',
      'COMPARE OBSERVATIONAL SUPPORT (original historical monthly mean vs field-year EOG):',
    ]
    for row in records:
        lines.append(f'{row["model"]}: n={row["n"]}, '
         f'yf={row["yf"]:.6f}, y0={row["y0"]:.6f}, k={row["k"]:.6f}, '
         f'R2={row["R2"]:.6f}, RMSE={row["RMSE"]:.6f}, '
         f'pixel CV={row["pixel_grouped_CV_RMSE"]:.6f}')
    lines += ['', 'Moth fit is a CANDIDATE: no changes to locked spider transfer.',
               'Review original moth PRC and R filename order and excluded Slabhouse QA,',
               'then lock independent moth transfer if accepted. This does not quantify',
               'actual behaviour or visible light under changing LED spectra.']
    (out/'moth_prc_fit_summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines),flush=True)
    print('Output:',out,flush=True)


if __name__=='__main__':
    main()
