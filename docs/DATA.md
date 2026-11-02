# Data and provenance

The Git repository contains source, compact canonical results, paper data,
and checksums. Large bases and proof inputs are companion Release assets.
Upload the files from the prepared companion-assets directory as Release
assets; do not add that directory to the source tree.

`data/archive_manifest.json` binds every asset and member by SHA-256, path,
and byte length. Members unpack relative to an explicitly chosen experiment
work directory. `data/provenance.json` records archived-source and public-copy
hashes. Public JSON references were made relative, and dependent file digests
were updated consistently. Self-protected preparation payload digests were
recomputed when necessary. Exact basis entries, rational norms, witnesses,
enumeration results, and historical native source/binary identities were
preserved. These transformations are packaging changes, not new searches.

Completed SQLite frontiers are read-only backups of queues with only completed
or split jobs. Their binary bytes may differ from the original database pages;
the public member hash binds the backup. Internal numerical job payloads and
completion records are retained. Running or partial queues are not included.

| Asset group | Contents |
| --- | --- |
| `inputs` | Seeded instances, original `B` bases, reduction transforms, original minima inputs |
| `projection` | Full and jet geometry, embeddings, saved minima and checked candidates |
| `bounds` | Original-lattice interval witnesses and exact lower-form certificates |
| `optimization` | Proposed/tested changes of basis, selected bases, block certificates |
| `frontiers` | Completed archived search-frontier snapshots |
| `index` | The selected 21-case exact-certificate index |

The canonical case JSON files keep squared lengths as exact rational numbers.
Their base-2 length logarithms are presentation values. `value_kind` identifies
archived exact values, checked candidate bounds, heuristic predictions, and
rounded paper reports. The five unfinished full projected values remain
candidate upper bounds; they must not be substituted for certified lower bounds.

The revised-paper upper bounds in `data/published/revised_lambda2_upper_bounds.csv`
were transcribed from printed Table 5 on page 8 (PDF page 9). The CSV records
the source document hash and printed precision. The source PDF is not bundled.
A rounded printed upper bound is a report of the source calculation, not an
exact rational certificate. Exact local upper-witness squared norms have their
own sources in `results/cases/`.

Large derived Bareiss files can be regenerated with `rebuild_preparation.py`;
the rebuilt checksum must agree. Full arrays, coefficient bounds, trusted
native source hashes, and completed-search status are part of the proof
provenance. Hash agreement alone is an integrity check, not a proof that an
enumeration algorithm is correct.
