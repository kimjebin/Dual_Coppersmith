#!/usr/bin/env python3
"""Unbounded preparation lifecycle and archive-scope regressions (no Sage)."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import signal
import sys
import tarfile
import tempfile
import threading
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from second_minimum import portable_minima, prepare_grid


class PreparationGridTests(unittest.TestCase):
    def test_no_case_timeout_and_complete_ordered_dimensions(self):
        rows = [('acd', 2, 1, 10), ('hm', 3, 3, 2)]
        waits = []

        class Process:
            def __init__(self, command, **options):
                self.pid = len(waits) + 1000
                folder = Path(options['stdout'].name).parent
                (folder / 'bundle.json').write_text(json.dumps(dict(
                    preparation_complete=True, ambient_lll_seconds=1, quotient_lll_seconds=2,
                    preparation_seconds=3, ambient_input_provenance='original_public_basis',
                    legacy_native_reused=False)))

            def wait(self, **options):
                waits.append(options)
                return 0

        with tempfile.TemporaryDirectory(prefix='prepare-grid-test-') as temporary:
            with patch.object(prepare_grid, 'GRID', rows), \
                    patch.object(prepare_grid.shutil, 'which', return_value='fake-sage'), \
                    patch.object(prepare_grid.subprocess, 'Popen', Process), redirect_stdout(io.StringIO()):
                code = prepare_grid.main(['--out', temporary, '--jobs', '2'])
            self.assertEqual(code, 0)
            self.assertEqual(waits, [{}, {}])
            saved = json.loads((Path(temporary) / 'preparations.json').read_text())
            self.assertEqual([row['enumeration_dimension'] for row in saved], [65, 9])
            manifest = json.loads((Path(temporary) / 'manifest.json').read_text())
            self.assertIsNone(manifest['case_timeout_seconds'])
            self.assertIsNone(manifest['preparation_timeout_seconds'])

    def test_worker_exception_cleans_all_started_groups(self):
        rows = [('hm', 3, 3, ell) for ell in (2, 3, 4, 5)]
        processes = []
        mutex = threading.Lock()
        killed = []

        class Process:
            def __init__(self, command, **options):
                self.event = threading.Event()
                with mutex:
                    self.pid = len(processes) + 1000
                    processes.append(self)

            def wait(self, **options):
                if self.pid == 1000:
                    return 0  # Missing bundle triggers a worker exception.
                if not self.event.wait(timeout=10):
                    raise AssertionError('active case was not cleaned up')
                return 130

        def kill_group(pid, signum):
            killed.append((pid, signum))
            next(process for process in processes if process.pid == pid).event.set()

        with tempfile.TemporaryDirectory(prefix='prepare-abort-test-') as temporary:
            with patch.object(prepare_grid, 'GRID', rows), \
                    patch.object(prepare_grid.shutil, 'which', return_value='fake-sage'), \
                    patch.object(prepare_grid.subprocess, 'Popen', Process), \
                    patch.object(prepare_grid.os, 'killpg', kill_group), redirect_stdout(io.StringIO()):
                with self.assertRaises(FileNotFoundError):
                    prepare_grid.main(['--out', temporary, '--jobs', '2'])
            for process in processes[1:]:
                self.assertIn((process.pid, signal.SIGINT), killed)
                self.assertIn((process.pid, signal.SIGKILL), killed)


class PortableExportTests(unittest.TestCase):
    def test_export_excludes_host_binaries_and_unrelated_root_files(self):
        with tempfile.TemporaryDirectory(prefix='portable-export-test-') as temporary:
            root = Path(temporary)
            grid = root / 'grid'
            case = grid / 'hm-t3-m3-ell2-seed0'
            (case / '.native').mkdir(parents=True)
            for filename in ('bundle.json', 'input.json', 'exact_enum-old'):
                (case / filename).write_text('{}')
            (case / '.native' / 'worker').write_text('host binary')
            (grid / 'unrelated.txt').write_text('not part of this export')
            (grid / 'manifest.json').write_text('{}')
            archive = root / 'bundle.tar.gz'
            with redirect_stdout(io.StringIO()):
                self.assertEqual(portable_minima.main(['export', str(grid), '--archive', str(archive)]), 0)
            with tarfile.open(archive) as incoming:
                names = incoming.getnames()
            self.assertEqual(sorted(names), sorted([
                'grid/hm-t3-m3-ell2-seed0/bundle.json',
                'grid/hm-t3-m3-ell2-seed0/input.json', 'grid/manifest.json']))

    def test_export_rejects_partial_grid_or_archive_inside_source(self):
        with tempfile.TemporaryDirectory(prefix='portable-export-input-test-') as temporary:
            root = Path(temporary)
            complete, incomplete = root / 'complete', root / 'incomplete'
            complete.mkdir()
            incomplete.mkdir()
            (complete / 'bundle.json').write_text('{}')
            (complete / 'input.json').write_text('{}')
            (incomplete / 'input.json').write_text('{}')
            with self.assertRaises(SystemExit):
                portable_minima.main(['export', str(root), '--archive', str(root / 'inside.tar.gz')])
            with self.assertRaises(SystemExit):
                portable_minima.main(['export', str(root), '--archive', str(root.parent / 'outside.tar.gz')])


if __name__ == '__main__':
    unittest.main()
