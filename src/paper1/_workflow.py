"""Explicit adapters around legacy implementations; no implicit YAML model changes."""
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/paper1"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance(script, args):
    version = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True)
    return dict(created_utc=datetime.now(timezone.utc).isoformat(),
                code_version=version.stdout.strip() or None,
                command=[sys.executable, str(script), *args],
                source_script_sha256=digest(script),
                config_hashes={str(p.relative_to(ROOT)): digest(p)
                               for p in sorted(CONFIG.glob("*")) if p.is_file()},
                config_usage="Recorded templates; explicit forwarded CLI arguments control legacy code.")


def forward(script_name, args=None):
    args = list(sys.argv[1:] if args is None else args)
    script = ROOT / "src" / script_name
    if not script.is_file():
        raise SystemExit("Missing existing implementation: " + str(script))
    # Help must not write sidecars or require local inputs.
    if "--help" in args or "-h" in args:
        raise SystemExit(subprocess.run([sys.executable, str(script), *args]).returncode)
    if "--output" not in args:
        raise SystemExit("Supply an explicit fresh --output directory")
    index = args.index("--output")
    if index + 1 >= len(args):
        raise SystemExit("--output requires a path")
    output = Path(args[index + 1])
    sidecar = output.parent / (output.name + ".paper1-provenance.json")
    if output.exists() or sidecar.exists():
        raise SystemExit("Output or provenance sidecar already exists; choose a fresh run")
    record = provenance(script, args)
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    record["status"] = "started"
    sidecar.write_text(json.dumps(record, indent=2) + "\n")
    code = subprocess.run(record["command"], cwd=Path.cwd()).returncode
    record.update(returncode=code, status="completed" if code == 0 else "failed")
    sidecar.write_text(json.dumps(record, indent=2) + "\n")
    raise SystemExit(code)


def pending(stage, specification):
    if "--help" in sys.argv[1:] or "-h" in sys.argv[1:]:
        print(stage + ": pending implementation. See " + specification)
        raise SystemExit(0)
    print(stage + ": PENDING; no scientific result produced. See " + specification,
          file=sys.stderr)
    raise SystemExit(2)
