import filecmp
import unittest
from pathlib import Path

from bin.adjust_pseudofinder_output import main


class TestAdjustPseudofinderOutput(unittest.TestCase):
    def setUp(self):
        self.prokka_standard = (
            Path(__file__).resolve().parent
            / "data/pseudofinder_test_inputs/prokka_standard.gff"
        )
        self.prokka_compliant = (
            Path(__file__).resolve().parent
            / "data/pseudofinder_test_inputs/prokka_compliant.gff"
        )
        self.pseudofinder_output = (
            Path(__file__).resolve().parent
            / "data/pseudofinder_test_inputs/pseudofinder_result.gff"
        )
        self.expected_output = (
            Path(__file__).resolve().parent
            / "data/pseudofinder_test_inputs/expected_processed_result.gff"
        )

    def test_main(self):
        test_output = "test_pseudofinder_processed_file.gff"
        main(
            self.pseudofinder_output,
            self.prokka_standard,
            self.prokka_compliant,
            test_output,
        )
        self.assertTrue(
            filecmp.cmp(test_output, self.expected_output, shallow=False),
            "Generated file does not match the expected file.",
        )
        # os.remove(test_output)
