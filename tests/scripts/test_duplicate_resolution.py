#!/usr/bin/env python3

# Copyright 2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
