import tempfile
import unittest
import csv
import shutil
from pathlib import Path

from reproduce import COMPARISONS, SOURCE_FILES, reproduce, mean, sem, welch_t_test, paired_t_test
from supplementary import statistics


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
            self.assertTrue((output / "supplementary_statistics.csv").is_file())

    def test_supplementary_p_values_and_scope(self):
        rows = statistics(ROOT / "data", mean, sem, welch_t_test, paired_t_test)
        self.assertEqual(len(rows), 28)
        self.assertEqual(sum(row["test"] == "Paired t test" for row in rows), 4)
        # Values independently checked against the final manuscript and figure exports.
        reported = {"S1A": .1963, "S1B": .0269, "S1C": .0353, "S1D": .0384, "S5A": .4959, "S5B": .0338, "S5C": .8068}
        for row in rows:
            if row["metric"] == "Total feeding time":
                self.assertEqual(round(row["p_value"], 4), reported[row["figure"]])
            else:
                self.assertGreater(row["p_value"], .05)

    def test_missing_paired_cell_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / "data"
            shutil.copytree(ROOT / "data" / "supplementary", copied / "supplementary")
            path = copied / "supplementary" / "Supplementary Figure S5" / "DAT-Cre_IBI.csv"
            with path.open(newline="") as handle:
                rows = list(csv.reader(handle))
            rows[2][1] = ""
            with path.open("w", newline="") as handle:
                csv.writer(handle).writerows(rows)
            with self.assertRaisesRegex(ValueError, "Missing or malformed"):
                statistics(copied, mean, sem, welch_t_test, paired_t_test)

    def test_changed_pair_order_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / "data"
            shutil.copytree(ROOT / "data" / "supplementary", copied / "supplementary")
            path = copied / "supplementary" / "mouse_values.csv"
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                fields = reader.fieldnames
                rows = list(reader)
            target = [r for r in rows if r["figure"] == "S5C" and r["group"] == "DAT-Cre CNO"]
            target[0]["pair_id"], target[1]["pair_id"] = target[1]["pair_id"], target[0]["pair_id"]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "pair IDs"):
                statistics(copied, mean, sem, welch_t_test, paired_t_test)


if __name__ == "__main__":
    unittest.main()
