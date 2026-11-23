"""Dispatch one explicit Paper 1 stage; dry-run unless --execute is supplied."""
import argparse
import shlex
import subprocess
import sys
from pathlib import Path

STAGES = {
    "audit": "audit_sites.py",
    "natural": "run_natural_patterns.py",
    "reconstruct-prc": "reconstruct_prc.py",
    "prc": "fit_pollution_prc.py",
    "luminance": "compare_equal_luminance.py",
    "conversion": "analyse_conversion_case.py",
    "fig01": "figures/fig01_natural_pattern.py",
    "fig02": "figures/fig02_pollution_patterns.py",
    "fig03": "figures/fig03_equal_prc.py",
    "fig04": "figures/fig04_equal_luminance.py",
    "fig05": "figures/fig05_conversion_case.py",
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", choices=STAGES, required=True)
    p.add_argument("--execute", action="store_true")
    p.add_argument("arguments", nargs=argparse.REMAINDER,
                   help="Legacy stage arguments after --")
    a = p.parse_args()
    args = a.arguments[1:] if a.arguments[:1] == ["--"] else a.arguments
    command = [sys.executable, str(Path(__file__).parent / STAGES[a.stage]), *args]
    print(shlex.join(command), flush=True)
    if a.execute:
        raise SystemExit(subprocess.run(command).returncode)
    print("Dry-run only; add --execute before -- to run this stage.")


if __name__ == "__main__":
    main()
