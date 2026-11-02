"""Locate the public baseline package without a machine-specific installation."""
from pathlib import Path
import os,sys
configured=os.environ.get("DUAL_COPPERSMITH_ROOT")
if configured:
    ROOT=Path(configured).resolve()
else:
    ROOT=next((root for root in Path(__file__).resolve().parents
               if (root/"shared/exact_minima.py").is_file()),None)
    if ROOT is None:raise RuntimeError("Cannot locate the public artifact root")
sys.path.insert(0,str(ROOT))

from lattices import *
HERE=Path(__file__).resolve().parent
