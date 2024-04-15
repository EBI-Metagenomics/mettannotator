#!/bin/env python3

import copy
import os
import re
import sys
import unittest

sys.path.append(os.path.sep.join(os.getcwd().split(os.path.sep)[:-2]))

from postprocessing.duplicate_resolution import main

DUMMY_REFERENCE = os.path.join(
    os.path.sep.join(os.getcwd().split(os.path.sep)[:-1]),
    "test_inputs",
    "dummy_reference.gff",
)
TEST_INPUT = os.path.join(
    os.path.sep.join(os.getcwd().split(os.path.sep)[:-1]),
    "test_inputs",
    "deduplication_script_test_input.gff",
)


class TestDuplicateResolution(unittest.TestCase):
    def test_main(self):
        result = main(
            DUMMY_REFERENCE, TEST_INPUT, "deduplication_test_outfile.gff", "uniformis"
        )
        self.assertEqual(result, {"replaced": 3, "unable_to_decide": 5})
