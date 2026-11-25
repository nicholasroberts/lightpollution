"""Copy a completed figure run into a fresh, checksum-manifested archive."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def archive_figure(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if not source.is_dir():
        raise ValueError('Figure output directory does not exist: '+str(source))
    if destination.exists():
        raise ValueError('Archive already exists; choose a new version: '+str(destination))
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError('Source and archive must be separate directory trees')
    files = sorted(p for p in source.rglob('*') if p.is_file())
    if not any(p.suffix.lower() in ('.png', '.pdf', '.svg') for p in files):
        raise ValueError('No figure artwork to archive')
    before = {str(p.relative_to(source)): digest(p) for p in files}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.mkdir(exist_ok=False)
    try:
        for path in files:
            target = destination/path.relative_to(source)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
        after = {name: digest(destination/name) for name in before}
        if before != after or any(digest(source/name) != sha for name, sha in before.items()):
            raise ValueError('Source changed or copy verification failed')
        (destination/'archive_manifest.json').write_text(json.dumps(dict(
            source_directory=str(source), files_sha256=after,
            archive_policy='Versioned snapshot; existing destinations are never overwritten. Retain source arrays at the paths and hashes in figure provenance.'
        ), indent=2)+'\n')
    except Exception:
        shutil.rmtree(destination)
        raise
    return destination


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--destination', type=Path, required=True)
    a = p.parse_args()
    try:
        print('Archived', archive_figure(a.source, a.destination))
    except (OSError, ValueError) as exc:
        p.error(str(exc))


if __name__ == '__main__':
    main()
