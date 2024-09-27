#!/bin/env python3

import os

from postprocessing.duplicate_resolution import main


def test_main(shared_datadir, tmp_path):
    result = main(
        (shared_datadir / "dummy_reference.gff"),
        (shared_datadir / "deduplication_script_test_input.gff"),
        (tmp_path / "deduplication_test_outfile.gff"),
        "uniformis",
    )
    os.remove("uniformis_replacements.txt")
    os.remove("uniformis_stats.txt")
    assert result == {"replaced": 6, "unable_to_decide": 5}
