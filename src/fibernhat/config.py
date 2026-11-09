"""Config = config/defaults.yaml <- config/site.yaml (gitignored) <- --config FILE.

Env FIBERNHAT_SITE_CONFIG points to an alternative site.yaml.
"""
from __future__ import annotations

import os
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]


def _merge(a: dict, b: dict) -> dict:
    out = dict(a)
    for k, v in (b or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_config(extra: str = None) -> dict:
    cfg = {"label_by": "auto", "threshold_estimate_reads": 20000, "filters": {}, "tools": {}}
    for p in (REPO_ROOT / "config" / "defaults.yaml",
              Path(os.environ.get("FIBERNHAT_SITE_CONFIG", REPO_ROOT / "config" / "site.yaml")),
              Path(extra) if extra else None):
        if p and p.exists():
            with open(p) as fh:
                cfg = _merge(cfg, yaml.safe_load(fh) or {})
    mk = cfg.get("tools", {}).get("modkit")
    if mk and not os.environ.get("FIBERNHAT_MODKIT"):
        os.environ["FIBERNHAT_MODKIT"] = str(mk)
    return cfg
