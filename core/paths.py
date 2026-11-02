"""Portable repository discovery for public commands and isolated workspaces."""
from pathlib import Path
import os

def repository_root():
    configured=os.environ.get("DUAL_COPPERSMITH_ROOT")
    if configured:
        root=Path(configured).resolve()
        if (root/"baseline/shared/exact_minima.py").is_file(): return root
        raise RuntimeError("DUAL_COPPERSMITH_ROOT is not a valid artifact root")
    for root in Path(__file__).resolve().parents:
        if (root/"baseline/shared/exact_minima.py").is_file(): return root
    raise RuntimeError("Cannot locate the artifact root")

ROOT=repository_root()
