"""Compatibility entrypoint for actual exact HM lambda2 with resumable workers."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from second_minimum.run_exact import main

options = [value for value in sys.argv[1:] if value != "--"]
if "--out" not in options and not any(value.startswith("--out=") for value in options):
    options += ["--out", os.path.join(ROOT, "second_minimum", "results", "exact-minima")]
raise SystemExit(main(["hm"] + options))
