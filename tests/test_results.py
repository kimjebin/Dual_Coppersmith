"""Guard numerical semantics and the archived/candidate distinction."""
from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from check import check_results,check_source

class ResultsTests(unittest.TestCase):
    def test_exact_rational_norms_intervals_and_table_regeneration(self):
        self.assertEqual(check_results()['unfinished'],5)
    def test_frozen_native_sources_and_python_syntax(self):
        self.assertEqual(check_source()['frozen_native_files'],4)

if __name__=='__main__':unittest.main()
