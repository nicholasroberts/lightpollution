#!/usr/bin/env python3
"""Madrid: UCM-local VIIRS versus full Robles et al. AstMon zenith observations.

Requires existing original Madrid municipality 2015-24 VNL V1 monthly CSV,
original 75N060W VNL V1 radiance/coverage GeoTIFFs, and optionally the
open-licensed Zenodo AstMon_evol.csv (doi:10.5281/zenodo.4633001).
No frozen annual, snow, VIIRS input rasters or PRC arrays are modified.

  python src/madrid_ucm_local_ground_bridge.py --self-test
  python src/madrid_ucm_local_ground_bridge.py --viirs
  python src/madrid_ucm_local_ground_bridge.py --download-ground --inspect-ground
  python src/madrid_ucm_local_ground_bridge.py --inspect-ground
  python src/madrid_ucm_local_ground_bridge.py --prepare-ground --ground-map CONFIG.json
  python src/madrid_ucm_local_ground_bridge.py --plot

The ground CSV's unknown original field/flag schema must be inspected before
specifying a CONFIG.json. In particular, never assume an integer flag means
'no moon, no clouds'. The mapping must explicitly identify accepted flag(s).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path
from tempfile import TemporaryDirectory

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Point

from local_viirs_five_case_analysis import ROOT, build_inventory, grid_month
from plot_madrid_viirs_robles_2021 import MONTHLY_CACHE

OUT = ROOT / 'data/processed/local_lighting_interventions/madrid_ground_viirs/local_ucm'
FIG = ROOT / 'figures/local_lighting_interventions/madrid_ground_viirs/local_ucm'
RAW = ROOT / 'data/raw/ground_sky_brightness/robles2021_zenodo4633001'
RAW_CSV = RAW / 'AstMon_evol.csv'
ZENODO = 'https://zenodo.org/records/4633001/files/AstMon_evol.csv?download=1'
ZENODO_MD5 = '634f098973251bf41f7d8a12b0a8ab45'
SITE = (-3 - 43/60 - 34/3600, 40 + 27/60 + 4/3600)  # UCM IAU I86
RADII_KM = (3, 10)  # Pre-specified, not selected by VIIRS outcome.
YM = pd.period_range('2015-01', '2019-12', freq='M').astype(str).tolist()
VIIRS_OUT = OUT / 'madrid_ucm_3km_10km_2015_2019_monthly.csv'
PAIR_OUT = OUT / 'madrid_ucm_fixed_common_pixels_2015_vs_2019.csv'
GROUND_MONTHLY_OUT = OUT / 'madrid_ucm_astmon_clear_after_midnight_2010_2020_monthly.csv'
GROUND_ANNUAL_OUT = OUT / 'madrid_ucm_astmon_clear_after_midnight_2010_2020_annual.csv'


def circle(lon, lat, radius_km):
    point = gpd.GeoSeries([Point(lon, lat)], crs='EPSG:4326').to_crs('EPSG:3035')
    return gpd.GeoSeries([point.iloc[0].buffer(radius_km * 1000)],
                         crs='EPSG:3035').to_crs('EPSG:4326').iloc[0]


def download_raw():
    RAW.mkdir(parents=True, exist_ok=True)
    if RAW_CSV.exists() and RAW_CSV.stat().st_size > 40_000_000:
        print('Existing AstMon CSV found; validating published MD5...')
    else:
        tmp = RAW_CSV.with_suffix('.csv.partial')
        request = urllib.request.Request(ZENODO, headers={'User-Agent':
            'Bristol lightpollution independent-validation/1.0'})
        try:
            with urllib.request.urlopen(request, timeout=120) as src, tmp.open('wb') as dest:
                while chunk := src.read(1024 * 1024):
                    dest.write(chunk)
        except Exception:
            tmp.unlink(missing_ok=True)
            raise
        tmp.replace(RAW_CSV)
    md5 = hashlib.md5()
    with RAW_CSV.open('rb') as src:
        while chunk := src.read(1024 * 1024):
            md5.update(chunk)
    if md5.hexdigest().lower() != ZENODO_MD5:
        raise ValueError(f'Zenodo published AstMon_evol.csv MD5 mismatch: '
                         f'{md5.hexdigest()} != {ZENODO_MD5}. '
                         'Do not analyse altered or partial ground data.')
    print('Verified Zenodo 4633001 AstMon CSV, MD5', md5.hexdigest())


def read_original(*, nrows=None):
    if not RAW_CSV.exists():
        raise FileNotFoundError(f'Download with --download-ground: {RAW_CSV}')
    with RAW_CSV.open('r', encoding='utf-8-sig', errors='replace') as f:
        sample = ''.join([f.readline() for _ in range(12)])
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=',;\t|')
        sep = dialect.delimiter
    except csv.Error:
        sep = ','
    frame = pd.read_csv(RAW_CSV, sep=sep, nrows=nrows, low_memory=False,
                        encoding='utf-8-sig')
    return frame


def inspect_ground():
    frame = read_original(nrows=200)
    print('Zenodo AstMon raw column names:', list(frame.columns))
    print('First 6 rows:\n', frame.head(6).to_string(index=False))
    for column in frame.columns:
        name = str(column).lower()
        if any(x in name for x in ('flag', 'cloud', 'moon', 'qual', 'clear', 'estado', 'status')):
            print(f'Possible quality field {column!r}:',
                  frame[column].astype(str).value_counts(dropna=False).head(15).to_dict())
    print('Full file must be analysed only with a reviewed mapping of '
          'date, B/V/R magnitude and *clear/moonless* flag(s). '
          'No automatic quality-code interpretation is permitted.')


def parse_datetime(values, spec):
    fmt = spec.get('time_format')
    if fmt == 'julian':
        jd = pd.to_numeric(values, errors='raise')
        return pd.to_datetime(jd, origin='julian', unit='D', utc=True)
    if fmt == 'iso_utc':
        return pd.to_datetime(values, errors='raise', utc=True, format='mixed')
    if fmt == 'iso_madrid_local':
        naive = pd.to_datetime(values, errors='raise', format='mixed')
        if getattr(naive.dt, 'tz', None) is not None:
            raise ValueError('Date has timezone; use iso_utc, not iso_madrid_local')
        return naive.dt.tz_localize(
            'Europe/Madrid', ambiguous='NaT', nonexistent='NaT').dt.tz_convert('UTC')
    raise ValueError('Mapping time_format must be julian, iso_utc or iso_madrid_local')


def prepare_ground(frame, spec, *, min_nights=3):
    """Parse only *explicitly approved* quality flags, then robust nightly medians."""
    if not spec.get('quality_filters'):
        raise ValueError('Mapping MUST include a nonempty quality_filters array '
                         'with documented clear/no-moon values. See --inspect-ground.')
    timecol = spec['time_col']
    if timecol not in frame:
        raise ValueError(f'Ground time column {timecol!r} missing')
    selected = pd.Series(True, index=frame.index)
    filter_counts = {}
    for rule in spec['quality_filters']:
        col, accepted = rule['column'], rule['accepted_values']
        if col not in frame or not accepted:
            raise ValueError(f'Unusable quality flag mapping {rule}')
        good = frame[col].astype(str).str.strip().str.casefold().isin(
            [str(v).strip().casefold() for v in accepted])
        selected &= good
        filter_counts[col] = {'accepted': list(accepted),
                              'matching_all_file': int(good.sum())}
    if selected.sum() < 100:
        raise ValueError('Too few accepted clear/no-moon observations; '
                         'check raw quality mappings, do not analyse unfiltered data')
    x = frame.loc[selected].copy()
    dt = parse_datetime(x[timecol], spec)
    x['madrid_datetime'] = dt.dt.tz_convert('Europe/Madrid')
    x = x.dropna(subset=['madrid_datetime'])
    x = x.loc[(x.madrid_datetime.dt.year.between(2010, 2020)) &
              (x.madrid_datetime.dt.hour >= 0) &
              (x.madrid_datetime.dt.hour < 5)]
    # 00:00-04:59 local; streetlights typically ramp back up ~04:45;
    # restrict further to strictly before 04:30.
    minute = 60*x.madrid_datetime.dt.hour + x.madrid_datetime.dt.minute
    x = x.loc[minute < 270].copy()
    if x.empty:
        raise ValueError('No observations after quality, midnight and time filters')
    x['night'] = x.madrid_datetime.dt.date
    x['year'] = x.madrid_datetime.dt.year
    x['month'] = x.madrid_datetime.dt.month
    records = []
    for band in ('B','V','R'):
        col = spec['bands'][band]
        if col not in x:
            raise ValueError(f'Ground magnitude column {col!r} missing')
        y = x[['year','month','night',col]].copy()
        y['mag'] = pd.to_numeric(y[col], errors='coerce')
        y = y.loc[y.mag.between(14, 24)]
        if y.empty:
            continue
        nights = (y.groupby(['year','month','night'], as_index=False)
                  .agg(night_median_mag=('mag','median'), n_obs=('mag','size')))
        months = (nights.groupby(['year','month'], as_index=False)
                  .agg(monthly_median_mag=('night_median_mag','median'),
                       night_count=('night','nunique'),
                       observations=('n_obs','sum')))
        months['band'] = band
        months.loc[months.night_count < min_nights, 'monthly_median_mag'] = np.nan
        months['month'] = months.month.astype(int)
        records.append(months)
    if len(records) != 3:
        raise ValueError('Missing one or more AstMon B/V/R bands after filters')
    monthly = pd.concat(records, ignore_index=True).sort_values(['year','month','band'])
    # Pooled per-observation median (paper used a different KDE/median method),
    # and separate equal-night-weight median to diagnose sampling effects.
    pooled = []
    for band in ('B','V','R'):
        y = x[['year',spec['bands'][band]]].copy()
        y['mag'] = pd.to_numeric(y[spec['bands'][band]], errors='coerce')
        y = y[y.mag.between(14, 24)]
        for year, group in y.groupby('year'):
            pooled.append({'year':int(year), 'band':band,
                           'pooled_obs_median_mag':float(group.mag.median()),
                           'n_observations':len(group)})
    annual = pd.DataFrame(pooled).sort_values(['year','band'])
    diagnostics = {'total_rows':len(frame), 'quality_rows':int(selected.sum()),
                   'after_midnight_rows':len(x), 'filters':filter_counts,
                   'ground_map':spec, 'source_md5':ZENODO_MD5,
                   'monthly_min_nights':min_nights,
                   'clock':'Europe/Madrid, local 00:00-04:29',
                   'notes':'The author used manual cloud/moon flags and KDE; '
                           'our ground aggregation is nightly and monthly MEDIANS, '
                           'not a reproduction of their Table 6 estimator.'}
    return monthly, annual, diagnostics


def viirs_near_ucm(index, *, years=range(2015, 2020), radii=RADII_KM):
    """Fixed centre-cell circular supports and strict paired-month pixel checks."""
    rowset, grids = [], {}
    for radius in radii:
        circle_poly = circle(*SITE, radius)
        key = f'UCM_{radius}km'
        shape_signature = None
        grids[key] = {}
        for year in years:
            for m in range(1,13):
                ym = f'{year}-{m:02d}'
                if ym not in index:
                    raise FileNotFoundError(f'EOG VCMSLCFG 75N060W missing {ym}')
                grid = grid_month(circle_poly, *index[ym])
                signature = (grid['crs'], grid['transform'],grid['rad'].shape)
                if shape_signature is None:
                    shape_signature = signature
                elif shape_signature != signature:
                    raise ValueError(f'{key} changed pixel grid {ym}')
                vals, cov = grid['rad'],grid['cf']
                valid = np.isfinite(vals)
                for threshold in (1, 3, 5):
                    mask = valid & (cov >= threshold)
                    rowset.append({
                        'support':key,'radius_km':radius,'year':year,
                        'month':m,'ym':ym,'cf_threshold':threshold,
                        'n_valid_pixels':int(mask.sum()),
                        'n_possible_pixels':int(grid['inside'].sum()),
                        'mean_radiance':(float(vals[mask].mean())
                                         if mask.any() else np.nan),
                        'median_cf_cvg':(float(np.median(cov[mask]))
                                         if mask.any() else np.nan),
                    })
                grids[key][ym] = (vals,cov)
            print(f'{key}: {year} 12/12 VNL V1 months', flush=True)
    paired = []
    for support, contents in grids.items():
        for m in range(1, 13):
            before = contents[f'2015-{m:02d}']
            after  = contents[f'2019-{m:02d}']
            for threshold in (1,3,5):
                mask = (np.isfinite(before[0]) & np.isfinite(after[0]) &
                        (before[1] >= threshold) & (after[1] >= threshold))
                n = int(mask.sum())
                av = float(before[0][mask].mean()) if n else np.nan
                bv = float(after[0][mask].mean()) if n else np.nan
                paired.append({'support':support, 'calendar_month':m,
                               'cf_threshold':threshold, 'paired_pixels':n,
                               'mean_2015':av, 'mean_2019':bv,
                               'delta':bv-av,
                               'percentage_delta':(100*(bv/av-1)
                                                   if n and av>0 else np.nan)})
    return pd.DataFrame(rowset),pd.DataFrame(paired)


def summary_2015_19(local, city, ground=None):
    result=[]
    for support, sub in local[local.cf_threshold.eq(1)].groupby('support'):
        means = sub.groupby('year').mean_radiance.mean()
        result.append({'series':support,'measurement':'upward VIIRS VNL V1',
                       'pct_2015_2019':100*(means.loc[2019]/means.loc[2015]-1),
                       'units':'within-series %'})
    for y in (2015,2019):
        if int((city.year.eq(y)).sum()) != 12:
            raise ValueError('Missing complete city-month 2015/19 controls')
    c = city.groupby('year').municipality_mean_vnlv1.mean()
    result.append({'series':'Madrid municipality',
                   'measurement':'upward VIIRS VNL V1',
                   'pct_2015_2019':100*(c.loc[2019]/c.loc[2015]-1),
                   'units':'within-series %'})
    if ground is not None:
        for band, sub in ground.groupby('band'):
            # EQUAL-calendar-month median magnitude, no imputation.
            for y in (2015,2019):
                if sub.loc[sub.year.eq(y),'monthly_median_mag'].notna().sum()<6:
                    print(f'WARNING: <6 valid ground months in {y} {band}, '
                          'omitting corresponding % contrast')
                    break
            else:
                annual = sub.groupby('year').monthly_median_mag.median()
                pct = 100*(10**(-.4*(annual.loc[2019]-annual.loc[2015]))-1)
                result.append({'series':f'AstMon Johnson {band}',
                               'measurement':'measured zenith sky brightness',
                               'pct_2015_2019':pct,
                               'units':'within-band %, median monthly mags'})
    return pd.DataFrame(result)


def plot_bridge(local, city, paired, ground=None, *, folder=FIG):
    folder.mkdir(parents=True,exist_ok=True)
    summary = summary_2015_19(local,city,ground)
    fig,axes=plt.subplots(3 if ground is not None else 2,1,
                          figsize=(12,9 if ground is not None else 6.7),
                          constrained_layout=True)
    a,b=axes[:2]
    viirs = local[local.cf_threshold.eq(1)].copy()
    for support,sub in viirs.groupby('support'):
        annual = sub.groupby('year').mean_radiance.mean()
        a.plot(annual.index,100*annual/annual.loc[2015],'-o',label=support)
    c = city.groupby('year').municipality_mean_vnlv1.mean()
    a.plot(c.index,100*c/c.loc[2015],'-o',label='Full municipality',
           color='black',linewidth=1.7)
    a.set(xlim=(2014.8,2019.2),xticks=range(2015,2020),
          ylabel='VIIRS VNL V1 radiance\n2015 = 100 for each area')
    a.set_title('A  Upward VIIRS: matched locations but different spatial supports')
    a.legend(ncol=3,fontsize=8); a.grid(axis='y',alpha=.2)
    for support,sub in paired[paired.cf_threshold.eq(3)].groupby('support'):
        b.plot(sub.calendar_month,sub.percentage_delta,'-o',label=support)
    b.axhline(0,color='black',linewidth=.8)
    b.set(xticks=range(1,13),xlabel='Calendar month',
          ylabel='2015→2019 VIIRS change (%)',
          title='B  Same-month, strictly paired VALID pixels (cf_cvg ≥3)')
    b.legend(fontsize=8); b.grid(axis='y',alpha=.2)
    if ground is not None:
        ax=axes[2]
        for band,sub in ground.groupby('band'):
            yearly=sub.groupby('year').monthly_median_mag.median().dropna()
            if 2015 not in yearly:
                continue
            vals=100*10**(-.4*(yearly-yearly.loc[2015]))
            ax.plot(yearly.index,vals,'-o',label=f'AstMon {band}',markersize=3)
        ax.axhline(100,color='black',lw=.7,ls='--')
        ax.set(xlim=(2014.8,2019.2),xticks=range(2015,2020),
               ylabel='Ground zenith brightness\n2015 = 100 in each band',
               title='C  UCM clear/moonless after-midnight sky brightness '
                     '(≥3 nights/month)')
        ax.legend(ncol=3,fontsize=8);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Madrid UCM observatory: local upward radiance and measured sky brightness\n'
                 'Fixed VIIRS footprints (3, 10 km) are NOT the ground skyglow source footprint',
                 fontsize=12.8,weight='bold')
    stem='madrid_ucm_local_viirs_2015_2019'
    if ground is not None:
        stem += '_and_astmon'
    for ext in ('.png','.svg','.pdf'):
        fig.savefig(folder/(stem+ext),bbox_inches='tight',dpi=180)
    plt.close(fig)
    return summary,stem


def self_test():
    assert abs(SITE[0]+3.726111111)<1e-6 and abs(SITE[1]-40.451111111)<1e-6
    # Test ground schema/flag gate, after midnight and per-night aggregator.
    stamp = pd.date_range('2015-01-01 00:10:00',periods=25,freq='D',tz='Europe/Madrid')
    fake = pd.DataFrame({'timestamp':stamp.repeat(5).tz_convert('UTC').astype(str),
                         'B':[19.27]*125,'V':[18.07]*125,'R':[17.42]*125,
                         'quality':['clear_nomoon']*100+['cloudy']*25})
    spec={'time_col':'timestamp','time_format':'iso_utc',
          'bands':{'B':'B','V':'V','R':'R'},
          'quality_filters':[{'column':'quality',
                              'accepted_values':['clear_nomoon']}]}
    monthly,yearly,log=prepare_ground(fake,spec)
    assert len(monthly)==3 and all(monthly.night_count==20)
    assert all(monthly.monthly_median_mag.notna())
    assert log['quality_rows']==100
    # Actual PNG/PDF/SVG plotting pipeline with 2015--19 synthetic data.
    period = pd.period_range('2015-01', '2019-12', freq='M')
    synthetic_local = pd.DataFrame([
        {'support':f'UCM_{r}km','cf_threshold':t,
         'year':int(v.year),'month':int(v.month),
         'mean_radiance':float(15+r/2-v.year/1000)}
        for r in RADII_KM for t in (1,3,5) for v in period])
    synthetic_city = pd.DataFrame({
        'year':period.year, 'calendar_month':period.month,
        'municipality_mean_vnlv1':np.linspace(30.,28.,len(period))})
    synthetic_pairs = pd.DataFrame([
        {'support':f'UCM_{r}km','cf_threshold':3,
         'calendar_month':m,'percentage_delta':-2.5}
        for r in RADII_KM for m in range(1,13)])
    with TemporaryDirectory(prefix='madrid-ucm-test-') as tmp:
        result, stem = plot_bridge(synthetic_local, synthetic_city,
                                    synthetic_pairs,folder=Path(tmp))
        assert len(result)==3
        assert all((Path(tmp)/(stem+ext)).exists()
                   for ext in ('.png','.pdf','.svg'))
    try:
        prepare_ground(fake,{**spec,'quality_filters':[]})
    except ValueError: pass
    else: raise AssertionError('Must refuse unfiltered ground data')
    # The fixed-radius circle must cover the paper's site and have plausible area.
    for radius in RADII_KM:
        poly=circle(*SITE,radius)
        assert poly.covers(Point(*SITE))
        area=float(gpd.GeoSeries([poly],crs='EPSG:4326')
                   .to_crs('EPSG:3035').area.iloc[0]/1e6)
        assert abs(area/(np.pi*radius**2)-1)<.025
    # Exercise real fixed-calendar-month common-pixel pairing with mock VIIRS.
    original_grid = globals()['grid_month']
    try:
        def tiny_grid(polygon, radfile, covfile):
            value = 10 if str(radfile).startswith('2015') else 8
            return {'rad': np.array([[value, value], [value, np.nan]],
                                    dtype=float),
                    'cf': np.array([[3, 5], [1, np.nan]],dtype=float),
                    'inside': np.array([[True, True],[True, True]]),
                    'crs':'EPSG:4326', 'transform':(0,1,0,0,0,-1)}
        globals()['grid_month'] = tiny_grid
        fake_index = {f'{year}-{month:02d}':(f'{year}-{month:02d}', 'coverage')
                      for year in (2015,2019) for month in range(1,13)}
        local, paired = viirs_near_ucm(fake_index, years=(2015,2019),
                                       radii=(3,))
        assert len(local)==72 and len(paired)==36
        p = paired[paired.cf_threshold.eq(3)]
        assert (p.paired_pixels==2).all()
        assert np.allclose(p.percentage_delta, -20)
    finally:
        globals()['grid_month'] = original_grid
    print('PASS: exact paper observatory coordinates, km-scale WGS84 geometry, '
          'fail-closed quality gate, full nightly/monthly AstMon aggregation')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--viirs',action='store_true',
                    help='Compute UCM 3/10km local 2015–19 raster series')
    ap.add_argument('--download-ground',action='store_true')
    ap.add_argument('--inspect-ground',action='store_true')
    ap.add_argument('--prepare-ground',action='store_true')
    ap.add_argument('--ground-map',type=Path,
                    help='Reviewed explicit JSON mapping of original Zenodo '
                         'time, B/V/R and clear/moonless quality flag values')
    ap.add_argument('--plot',action='store_true')
    args=ap.parse_args()
    if args.self_test:
        self_test();return
    if args.download_ground:
        download_raw()
    if args.inspect_ground:
        inspect_ground()
    if args.viirs:
        index,missing=build_inventory()
        if missing:
            raise FileNotFoundError(f'Incomplete European VIIRS 75N060W pair set: '
                                    f'{missing[:8]}')
        local,paired=viirs_near_ucm(index)
        OUT.mkdir(parents=True,exist_ok=True)
        local.to_csv(VIIRS_OUT,index=False)
        paired.to_csv(PAIR_OUT,index=False)
        print('Local UCM VIIRS monthly:',VIIRS_OUT)
        print('Strict 2015/2019 paired pixels:',PAIR_OUT)
    if args.prepare_ground:
        if not args.ground_map or not args.ground_map.exists():
            raise ValueError('Review --inspect-ground output then provide '
                             '--ground-map path with explicit quality flags')
        spec=json.loads(args.ground_map.read_text(encoding='utf-8'))
        ground,annual,log=prepare_ground(read_original(),spec)
        OUT.mkdir(parents=True,exist_ok=True)
        ground.to_csv(GROUND_MONTHLY_OUT,index=False)
        annual.to_csv(GROUND_ANNUAL_OUT,index=False)
        (OUT/'astmon_filter_audit.json').write_text(
            json.dumps(log,indent=2,default=str),encoding='utf-8')
        print('AstMon clear/moonless monthly:',GROUND_MONTHLY_OUT)
        print('Original flag selection audit:',OUT/'astmon_filter_audit.json')
    if args.plot or args.viirs or args.prepare_ground:
        if not VIIRS_OUT.exists() or not PAIR_OUT.exists():
            if args.viirs:
                raise ValueError('Requested VIIRS run did not generate data')
            print('First run --viirs to extract local 3/10 km supports');return
        if not MONTHLY_CACHE.exists():
            raise FileNotFoundError(f'Previously calculated Madrid municipality '
                                    f'CSV not found: {MONTHLY_CACHE}')
        local=pd.read_csv(VIIRS_OUT)
        paired=pd.read_csv(PAIR_OUT)
        city=pd.read_csv(MONTHLY_CACHE)
        ground=pd.read_csv(GROUND_MONTHLY_OUT) if GROUND_MONTHLY_OUT.exists() else None
        summary,stem=plot_bridge(local,city,paired,ground)
        summary.to_csv(OUT/'madrid_ucm_2015_2019_footprint_ground_comparison.csv',
                       index=False)
        print(summary.to_string(index=False))
        print('Figure:',FIG/(stem+'.png'))
    if not any((args.download_ground,args.inspect_ground,args.viirs,
                args.prepare_ground,args.plot)):
        ap.print_help()


if __name__=='__main__':
    main()