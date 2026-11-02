"""Compatibility entrypoint for actual exact ACD lambda2, including dense rows."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from second_minimum.run_exact import main

arguments = [value for value in sys.argv[1:] if value != "--"]
if len(arguments) < 3:
    raise SystemExit("Usage: sage acd_second_minimum.sage t m ell [--workers N --out DIR ...]")
t, m, ell = arguments[:3]
options = arguments[3:]
if "--out" not in options and not any(value.startswith("--out=") for value in options):
    options += ["--out", os.path.join(ROOT, "second_minimum", "results", "exact-minima")]
raise SystemExit(main(["acd", "--t", t, "--m", m, "--ell", ell] + options))
