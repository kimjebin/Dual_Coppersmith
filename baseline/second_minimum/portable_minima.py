#!/usr/bin/env python3
"""Verify portable case/grid bundles, or export an archive without host binaries."""

import argparse
from contextlib import closing, ExitStack
import json
from pathlib import Path
import sys
import tarfile
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.exact_minima import _lock, verify_minima_bundle


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=('verify', 'export'))
    parser.add_argument('directory', type=Path)
    parser.add_argument('--archive', type=Path)
    args = parser.parse_args(argv)
    directory = args.directory.resolve()
    cases = ([directory] if (directory / 'bundle.json').exists() else
             sorted(path.parent for path in directory.glob('*/bundle.json')))
    if not cases:
        parser.error('no completed portable bundles found')
    if args.operation == 'export' and args.archive is None:
        parser.error('export requires --archive (new .tar.gz file)')
    if args.operation == 'export':
        if args.archive.resolve().is_relative_to(directory):
            parser.error('archive must be outside the input directory')
        incomplete = [path.parent.name for path in directory.glob('*/input.json')
                      if not (path.parent / 'bundle.json').exists()]
        if incomplete:
            parser.error('incomplete cases cannot be exported: ' + ', '.join(incomplete))
    with ExitStack() as stack:
        for case in cases:
            stack.enter_context(_lock(case))
        if args.operation == 'verify':
            for case in cases:
                bundle = verify_minima_bundle(case)
                print(json.dumps(dict(case=case.name, verified=True,
                                      dimension=bundle['enumeration_dimension'],
                                      bundle_payload_sha256=bundle['payload_sha256'])), flush=True)
        else:
            # Locks keep committed frontiers stable. SQLite backup consolidates
            # WALs in the export only; runtime checkpoints are never modified.
            import sqlite3
            with tempfile.TemporaryDirectory(prefix='minima-export-') as temporary:
                with tarfile.open(args.archive.resolve(), 'x:gz', compresslevel=1) as archive:
                    paths = [path for case in cases for path in case.rglob('*')]
                    if cases != [directory]:
                        paths += [directory / name for name in
                                  ('manifest.json', 'preparations.json', 'preparations.csv', 'verification.json')
                                  if (directory / name).exists()]
                    for path in sorted(paths):
                        relative = path.relative_to(directory)
                        if (not path.is_file() or '.native' in relative.parts or
                                path.name == '.lock' or path.name.startswith(('exact_enum-', 'frontier-import-')) or
                                path.name.endswith(('.tmp', '-wal', '-shm'))):
                            continue
                        source = path
                        if path.name == 'frontier.sqlite3':
                            source = Path(temporary) / (path.parent.parent.name + '-frontier.sqlite3')
                            with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as incoming, \
                                    closing(sqlite3.connect(source)) as outgoing:
                                incoming.backup(outgoing)
                                outgoing.execute('PRAGMA journal_mode=DELETE')
                        archive.add(source, arcname=str(Path(directory.name) / relative), recursive=False)
            print(f'Exported {len(cases)} cases to {args.archive}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
