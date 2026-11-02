# Upstream software

Upstream author names are attribution required by the included licenses. They
do not identify the authors of this submission and have not been anonymized.

| Component | Source and license |
| --- | --- |
| Baseline small-root utilities | [crypto-attacks](https://github.com/jvdsn/crypto-attacks), MIT; notice in `baseline/shared/LICENSE.crypto-attacks` |
| Generated enumlib headers | MIT notices retained in both files in `experiments/engine/generated/` |
| fplll | [fplll](https://github.com/fplll/fplll), pinned commit and license in `manifest.json` and `fplll-5.5.0.LICENSE` |
| G6K (optional sieving) | [G6K](https://github.com/fplll/g6k), pinned commit and license in `manifest.json` and the accompanying G6K license file |

No upstream Git history, installation prefix, compiled library, or build cache
is bundled. The recorded local upstream checkouts had no tracked modifications.
The optional G6K build script checks out the pinned commit in a user-selected
work directory. Native enumeration links to the separately installed fplll,
GMP, and MPFR libraries.
