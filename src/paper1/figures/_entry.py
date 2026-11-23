"""Shared figure adapter entry points."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _workflow import forward, pending
