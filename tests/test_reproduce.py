import tempfile
import unittest
from pathlib import Path

from reproduce import COMPARISONS, SOURCE_FILES, reproduce


ROOT = Path(__file__).resolve().parents[1]


class ReproductionTests(unittest.TestCase):
    def test_exact_source_set(self):
        self.assertEqual(
            {path.name for path in (ROOT / "data").glob("*.csv")},
            set(SOURCE_FILES),
        )
    def test_all_manuscript_checks_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            result = reproduce(ROOT / "data", Path(directory) / "results")
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(int(result["regression"]["n"]), 170)
        self.assertEqual(round(result["regression"]["r_squared"], 4), 0.5819)
        self.assertLess(result["regression"]["p_value"], 0.0001)
        self.assertEqual(len(result["comparisons"]), len(COMPARISONS))
        self.assertTrue(all(row["status"] == "PASS" for row in result["comparisons"]))

    def test_outputs_are_created(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "results"
            reproduce(ROOT / "data", output)
            self.assertTrue((output / "manuscript_statistics.csv").is_file())
            self.assertTrue((output / "Figure_1D_regression.svg").is_file())
            self.assertTrue((output / "reproduction_report.md").is_file())
            self.assertTrue((output / "source_checksums.sha256").is_file())


if __name__ == "__main__":
    unittest.main()
