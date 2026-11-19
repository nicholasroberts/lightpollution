"""Register prepared RAW batches and build a local visual-review page.

Uses the existing single-site landmark registration without polarization fitting.
Review decisions are exported separately and never modify numerical locks.
"""
import argparse
import hashlib
import html
import json
import re
import subprocess
import sys
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_review(output, summary):
    esc = html.escape
    sections = []
    for row in summary['sites']:
        site = row['site']
        destination = output / 'runs' / site / site
        checks = row.get('exposures', [])
        table = ''.join('<tr><td>{}</td><td>{}</td><td>{}</td></tr>'.format(
            esc(t['filename']), 'PASS' if t['quality'].get('numerical_pass') else 'PENDING',
            esc(str(t['quality'].get('heldout_median_all_matches_px', '—')))) for t in checks)
        images = ''
        for name, label in [('registration_edges.png', 'Exposure alignment'),
                            ('registered_dop_and_north.png', 'Reference Moon and reconstructed DoLP')]:
            image = destination / name
            if image.exists():
                rel = image.relative_to(output).as_posix()
                images += f'<p>{label}: click to inspect at full resolution.</p><a href="{rel}" target="_blank"><img loading="lazy" src="{rel}" alt="{esc(site)} {label}"></a>'
        sections.append(f'''<section data-site="{esc(site)}"><h2>{esc(site)}</h2>
<p>{esc(row['status'])} · North: {esc(str(row.get('north_orientation_status', 'unresolved')))}</p>
<p>{esc(str(row.get('reason') or ''))}</p>
<label>Alignment <select class="alignment"><option>unreviewed</option><option>acceptable</option><option>needs work</option></select></label>
<label>Moon identity <select class="moon"><option>unreviewed</option><option>confirmed</option><option>uncertain</option><option>incorrect</option></select></label>
<label>Notes <input class="notes" placeholder="Exposure numbers / visible mismatch"></label>
<details><summary>Numerical checks ({len(checks)} exposures)</summary><table><tr><th>Exposure</th><th>Numerical</th><th>Held-out median (px)</th></tr>{table}</table></details>{images}</section>''')
    payload = json.dumps({'registration_summary_sha256': digest(output / 'registration_batch_summary.json'),
                          'status': 'visual_review_only_not_model_acceptance'})
    page = '''<!doctype html><meta charset="utf-8"><title>Green RAW alignment review</title>
<style>body{font:16px system-ui;max-width:1500px;margin:24px auto;padding:0 20px;background:#fafafa}header{position:sticky;top:0;background:#fff;padding:12px;z-index:2;border-bottom:1px solid #aaa}section{background:white;padding:20px;margin:20px 0;border:1px solid #ccc}img{width:100%;height:auto}label{display:inline-block;margin:8px}input{width:350px}td,th{padding:5px 18px;text-align:left}button{margin:6px;padding:8px}select,input{font:inherit}</style>
<header><b>Alignment review</b> <button id="previous">Previous site</button><button id="next">Next site</button><button id="export">Download review JSON</button><label>Restore review JSON <input id="restore" type="file" accept=".json"></label><span id="count"></span></header>
<h1>Green RAW registration</h1><p>Red = reference; cyan = aligned exposure; white = overlapping edges. Inspect buildings and skyline across the full circle, then stars where useful. Look for doubled edges, especially in the last analyzer group. Moving stars and the Moon need not overlap exactly across time.</p>
<p>Numerical PASS does not establish visual acceptance. Moon identity is reviewed separately. North and projection remain provisional; this page does not validate handedness or angular calibration. Download your decisions before closing; no decisions are saved automatically.</p>
''' + '\n'.join(sections) + '''<script>
const base=PAYLOAD;
const cards=Array.from(document.querySelectorAll('section'));
let index=0;
function go(delta){index=Math.max(0,Math.min(cards.length-1,index+delta));cards[index]?.scrollIntoView({block:'start'});}
document.querySelector('#previous').onclick=()=>go(-1);
document.querySelector('#next').onclick=()=>go(1);
function update(){document.querySelector('#count').textContent=cards.filter(c=>c.querySelector('.alignment').value!=='unreviewed').length+' / '+cards.length+' reviewed';}
document.addEventListener('change',update);update();
document.querySelector('#export').onclick=()=>{const sites=cards.map(c=>({site:c.dataset.site,alignment:c.querySelector('.alignment').value,moon_identity:c.querySelector('.moon').value,notes:c.querySelector('.notes').value}));const b=new Blob([JSON.stringify({...base,sites},null,2)],{type:'application/json'});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='alignment_visual_review.json';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);};
document.querySelector('#restore').onchange=async e=>{try{const r=JSON.parse(await e.target.files[0].text());if(r.registration_summary_sha256!==base.registration_summary_sha256)throw Error('Review belongs to a different registration summary');for(const c of cards){const v=r.sites.find(s=>s.site===c.dataset.site);if(v){c.querySelector('.alignment').value=v.alignment;c.querySelector('.moon').value=v.moon_identity;c.querySelector('.notes').value=v.notes||'';}}update();}catch(err){alert(err.message);}};
</script>'''.replace('PAYLOAD', payload)
    (output / 'review.html').write_text(page, encoding='utf-8')


def run_batch(a):
    source = json.loads(a.raw_batch.read_text())
    if source.get('status') != 'conversion_complete_not_registered':
        raise ValueError('Require a complete RAW conversion batch')
    rows = source['sites']
    if a.sites:
        missing = set(a.sites) - {r['site'] for r in rows}
        if missing:
            raise ValueError('Unknown sites: ' + ', '.join(sorted(missing)))
        rows = [r for r in rows if r['site'] in a.sites]
    if not rows or len({r['site'] for r in rows}) != len(rows):
        raise ValueError('Empty or duplicate site list')
    # Preflight every manifest before starting any expensive registration.
    for r in rows:
        if not re.fullmatch(r'[A-Za-z0-9_-]+', r['site']):
            raise ValueError('Unsafe site name')
        if r['status'] != 'prepared_not_registered' or r['returncode'] != 0:
            raise ValueError('Conversion failed: ' + r['site'])
        if digest(r['manifest_path']) != r['manifest_sha256']:
            raise ValueError('Prepared manifest changed: ' + r['site'])
    if not (a.model_run / 'hps/provenance.json').is_file():
        raise ValueError('Model provenance missing')
    a.output.mkdir(parents=True, exist_ok=False)
    (a.output / 'runs').mkdir()
    summary = dict(status='registration_in_progress', camera_channel=source['camera_channel'],
                   conversion_batch_sha256=digest(a.raw_batch), sites=[], visual_review_required=True,
                   sky_pattern_fit_performed=False, raw_inputs_modified=False)
    script = Path(__file__).with_name('lock_spider_site_registration.py')
    for index, r in enumerate(rows):
        site = r['site']
        row = dict(site=site, status='pending_inputs_or_registration')
        print(f'[{index+1}/{len(rows)}] {site}', flush=True)
        report = a.batch / site / 'registration/rotation_analysis.json'
        if not report.is_file():
            row['reason'] = 'Missing prior reference geometry/report: ' + str(report)
        else:
            destination = a.output / 'runs' / site
            command = [sys.executable, str(script), '--batch', str(a.batch), '--model-run', str(a.model_run),
                       '--sites', site, '--registration-model', 'fisheye-shared',
                       '--linear-raw-manifest', r['manifest_path'], '--output', str(destination)]
            if a.moon_coordinates:
                command += ['--moon-coordinates', str(a.moon_coordinates)]
            with (a.output / (site + '_registration.log')).open('w') as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
            row.update(returncode=result.returncode, log=site + '_registration.log')
            lock = destination / 'locked_site_references.json'
            if lock.exists():
                record = json.loads(lock.read_text())['sites'][site]
                row.update(status=record['status'], reason=record.get('reason'),
                           north_orientation_status=record.get('north_orientation_status'),
                           moon_localization=record.get('moon_localization'),
                           exposures=[dict(filename=t['filename'], quality=t['quality']) for t in record.get('transforms', [])],
                           lock_path=lock.relative_to(a.output).as_posix(), lock_sha256=digest(lock))
            if result.returncode:
                row.update(status='registration_process_failed', reason='See registration log')
        summary['sites'].append(row)
        (a.output / 'registration_batch_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        write_review(a.output, summary)
        print('  ' + row['status'] + ': ' + str(row.get('reason') or 'review overlays'), flush=True)
    summary['status'] = 'registration_attempts_complete_visual_review_pending'
    (a.output / 'registration_batch_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    write_review(a.output, summary)
    print('Open ' + str(a.output / 'review.html'), flush=True)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-batch', type=Path, required=True, help='Prepared conversion batch_summary.json')
    p.add_argument('--batch', type=Path, required=True, help='Original TIFF registration batch with reference geometry')
    p.add_argument('--model-run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sites', nargs='+')
    p.add_argument('--moon-coordinates', type=Path)
    a = p.parse_args()
    try:
        summary = run_batch(a)
    except (OSError, ValueError, KeyError) as exc:
        p.error(str(exc))
    return 0 if all(r['status'] == 'numerically_locked_awaiting_visual_review' for r in summary['sites']) else 1


if __name__ == '__main__':
    sys.exit(main())
