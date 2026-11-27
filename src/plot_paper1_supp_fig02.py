"""Render Supplementary Figure 2 from saved Oxford fits, without refitting."""
import argparse
import json
from pathlib import Path
import shutil
from run_paper1_oxford_natural_fits import SITES, digest, render

def produce(root, output):
    if output.exists():
        raise ValueError('Use a fresh version directory; existing figures are retained')
    sources = [root/'site_fit_summary.json', root/'site_fit_summary.csv']
    for site in SITES:
        base=root/site
        sources.extend(base/'fit'/name for name in (
            'solar_lunar_fit.json','fit_diagnostics.png',
            'ratio_profile.csv','separation_profile.csv'))
    for path in sources:
        if not path.is_file():
            raise ValueError('Missing saved result: '+str(path))
    # Figure 1 renderer verifies source-lock, north-data and model provenance,
    # recomputes the saved components and checks the saved mixture before plotting.
    render(root, output, supplemental=True)
    copied=[]
    for source in sources:
        dest=output/'results'/source.relative_to(root)
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,dest)
        copied.append(dict(source=str(source),copy=str(dest.relative_to(output)),
                           sha256=digest(dest)))
    report=dict(figure='Supplementary Figure 2',sites=list(SITES),
        fit_repeated=False,source_root=str(root),
        layout='Four rows: DoLP/AoP per site; measured/lunar/solar/combined columns',
        colormaps={'DoLP':'viridis linear 0–1','AoP':'romaO cyclic −90–90 degrees'},
        files=copied,script_sha256=digest(__file__))
    (output/'archive_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Saved Supplementary Figure 2:',output)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('data/processed/paper1_natural_fit_oxford_01'))
    p.add_argument('--output',type=Path,default=Path('figures/paper1/supplemental/supp_fig02/v01'))
    a=p.parse_args()
    try:produce(a.root,a.output)
    except (OSError,ValueError,KeyError) as exc:p.error(str(exc))

if __name__=='__main__':main()
