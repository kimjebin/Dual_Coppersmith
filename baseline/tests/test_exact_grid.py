#!/usr/bin/env python3
"""Whole-case timing and process-group deadline regressions."""

import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from second_minimum.run_exact_grid import Campaign, PAPER_CASES, run_case, stop_group, write_table
from second_minimum.run_exact_grid import atomic_json


class GridTests(unittest.TestCase):
    def test_timeout_stops_descendants_even_when_leader_exits(self):
        with tempfile.TemporaryDirectory(prefix="exact-grid-child-test-") as output:
            pidfile = Path(output) / "child.pid"
            child = "import signal,time;signal.signal(signal.SIGINT,signal.SIG_IGN);time.sleep(30)"
            parent = ("import subprocess,sys,time;from pathlib import Path;"
                      f"p=subprocess.Popen([sys.executable,'-c',{child!r}]);"
                      f"Path({str(pidfile)!r}).write_text(str(p.pid));time.sleep(30)")
            proc = subprocess.Popen([sys.executable, "-c", parent], start_new_session=True,
                                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            try:
                deadline = time.monotonic() + 5
                while not pidfile.exists() and time.monotonic() < deadline:
                    time.sleep(0.01)
                self.assertTrue(pidfile.exists())
                pid = int(pidfile.read_text())
                stop_group(proc, grace=0.1)
                self.assertIsNotNone(proc.poll())
                for _ in range(100):
                    stat = Path(f"/proc/{pid}/stat")
                    if not stat.exists() or stat.read_text().split()[2] == "Z":
                        break
                    time.sleep(0.01)
                else:
                    self.fail("Descendant remained running after deadline cleanup")
            finally:
                if proc.poll() is None:
                    stop_group(proc, grace=0.1)

    def test_run_case_reports_outer_deadline_without_claiming_exact(self):
        real_popen = subprocess.Popen
        def fake_sage(command, **options):
            return real_popen([sys.executable, "-c", "import time;time.sleep(30)"], **options)
        with tempfile.TemporaryDirectory(prefix="exact-grid-timeout-test-") as output:
            with patch("second_minimum.run_exact_grid.subprocess.Popen", side_effect=fake_sage):
                row = run_case(("hm", 3, 3, 2), output=Path(output), timeout=0.1,
                               workers=1, memory_gb=1)
            self.assertEqual(row["status"], "timeout")
            self.assertFalse(row["exact"])
            self.assertFalse(row["enumeration_entered"])
            self.assertIsNone(row["lattice_dimension"])
            self.assertGreaterEqual(row["elapsed_seconds"], 0.1)
            self.assertIsNotNone(row["returncode"])
            folder = Path(row["checkpoint"])
            self.assertEqual(json.loads((folder / "timing.json").read_text())["status"], "timeout")

    def test_campaign_abort_cleans_active_group_and_prevents_late_start(self):
        campaign = Campaign()
        proc = campaign.start([sys.executable, "-c", "import time;time.sleep(30)"],
                               start_new_session=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
        campaign.abort()
        self.assertIsNotNone(proc.poll())
        with self.assertRaises(InterruptedError):
            campaign.start([sys.executable, "-c", "pass"], start_new_session=True)
        campaign.release(proc)
        self.assertEqual(campaign.active, set())

    def test_metadata_write_failure_cleans_started_process(self):
        real_popen = subprocess.Popen
        processes = []
        def fake_sage(command, **options):
            proc = real_popen([sys.executable, "-c", "import time;time.sleep(30)"], **options)
            processes.append(proc)
            return proc
        with tempfile.TemporaryDirectory(prefix="exact-grid-write-test-") as output:
            with patch("second_minimum.run_exact_grid.subprocess.Popen", side_effect=fake_sage), \
                 patch("second_minimum.run_exact_grid.atomic_json", side_effect=OSError("disk full")):
                with self.assertRaises(OSError):
                    run_case(("hm", 3, 3, 2), output=Path(output), timeout=10,
                             workers=1, memory_gb=1)
        self.assertEqual(len(processes), 1)
        self.assertIsNotNone(processes[0].poll())

    def test_exit_observed_after_deadline_is_still_timeout(self):
        real_popen = subprocess.Popen
        def fake_sage(command, **options):
            return real_popen([sys.executable, "-c", "pass"], **options)
        def delayed_metadata(path, value):
            atomic_json(path, value)
            if path.name == "supervisor.json":
                time.sleep(0.2)
        with tempfile.TemporaryDirectory(prefix="exact-grid-late-exit-test-") as output:
            with patch("second_minimum.run_exact_grid.subprocess.Popen", side_effect=fake_sage), \
                 patch("second_minimum.run_exact_grid.atomic_json", side_effect=delayed_metadata):
                row = run_case(("hm", 3, 3, 2), output=Path(output), timeout=0.1,
                               workers=1, memory_gb=1)
            self.assertEqual(row["returncode"], 0)
            self.assertEqual(row["status"], "timeout")
            self.assertFalse(row["exact"])

    def test_table_preserves_grid_order_and_timeout_status(self):
        rows = [{"application": "HM", "t": 3, "m": 3, "ell": 9, "status": "timeout"},
                {"application": "ACD", "t": 2, "m": 1, "ell": 10, "status": "complete"}]
        with tempfile.TemporaryDirectory(prefix="exact-grid-table-test-") as output:
            ordered = write_table(Path(output), rows)
            self.assertEqual(ordered[0]["application"], "ACD")
            saved = json.loads((Path(output) / "timings.json").read_text())
            self.assertEqual(saved[1]["status"], "timeout")
        self.assertEqual(len(PAPER_CASES), 21)
        self.assertEqual(len(set(PAPER_CASES)), 21)


if __name__ == "__main__":
    unittest.main()
