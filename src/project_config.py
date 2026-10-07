"""Configuration helpers for the lightpollution pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "pipeline.yaml"


def load_config(path: Path | str = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load the YAML project configuration."""
    config_path = Path(path)
    if not config_path.is_absolute():
        config_path = PROJECT_ROOT / config_path

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise ValueError(f"Configuration is not a mapping: {config_path}")

    return config


def project_path(relative_path: str | Path) -> Path:
    """Resolve a path relative to the repository root."""
    return PROJECT_ROOT / Path(relative_path)


def configured_path(config: dict[str, Any], key: str) -> Path:
    """Resolve one entry from the config paths section."""
    try:
        relative = config["paths"][key]
    except KeyError as exc:
        raise KeyError(f"Missing paths.{key} in pipeline configuration") from exc

    return project_path(relative)
