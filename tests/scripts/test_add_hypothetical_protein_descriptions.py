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

import copy
import re

import pytest

from bin.add_hypothetical_protein_descriptions import (
    GeneCaller,
    clean_up_function,
    escape_reserved_characters,
    get_function,
    insert_product_source,
    keep_or_move_to_note,
    move_function_to_note,
)


@pytest.fixture
def setup_data():
    eggnog_annot = {"my_id": "eggnog function"}
    attributes_dict = {
        "ID": "my_id",
        "product": "hypothetical protein",
        "protein_id": "protein id",
    }
    ipr_info = {
        "my_id": {
            "Family": {
                "NCBIFam": {
                    "match": 0.80,
                    "ipr_desc": "NCBI Fam function",
                    "sig_desc": "sig desc",
                    "level": 1,
                },
                "Pfam": {
                    "match": 0.90,
                    "ipr_desc": "Pfam Fam function",
                    "sig_desc": "sig desc",
                    "level": 0,
                },
            },
            "Domain": {
                "CDD": {
                    "match": 1.00,
                    "ipr_desc": "CDD Domain function",
                    "sig_desc": "sig desc",
                    "level": 0,
                }
            },
        }
    }
    ipr_memberdb_only = {
        "my_id": {
            "no_type": {
                "NCBIFam": {
                    "match": 0.80,
                    "ipr_desc": "NCBI Fam function",
                    "sig_desc": "sig desc",
                    "level": 2,
                },
                "Pfam": {
                    "match": 0.90,
                    "ipr_desc": "Pfam Fam function",
                    "sig_desc": "sig desc",
                    "level": 0,
                },
            },
        }
    }
    return eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only


def test_get_function_ncbifam_priority(setup_data):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    result = get_function(
        "my_id",
        attributes_dict,
        eggnog_annot,
        ipr_info,
        ipr_memberdb_only,
        GeneCaller.PROKKA,
    )
    assert result == ("NCBI Fam function", "InterPro(NCBIFam)")


def test_get_function_level_priority(setup_data):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    ipr_test = copy.deepcopy(ipr_info)
    ipr_test["my_id"]["Family"]["Pfam"].update({"level": 1})
    result = get_function(
        "my_id",
        attributes_dict,
        eggnog_annot,
        ipr_test,
        ipr_memberdb_only,
        GeneCaller.PROKKA,
    )
    assert result == ("Pfam Fam function", "InterPro(Pfam)")


def test_get_function_level_eggnog(setup_data):
    eggnog_annot, attributes_dict, _, _ = setup_data
    result = get_function(
        "my_id", attributes_dict, eggnog_annot, {}, {}, GeneCaller.PROKKA
    )
    assert result == ("eggnog function", "eggNOG")


def test_get_function_level_pfam_overtake(setup_data):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    ipr_test = copy.deepcopy(ipr_info)
    del ipr_test["my_id"]["Family"]["NCBIFam"]
    ipr_test["my_id"]["Family"].update({"CDD": ipr_test["my_id"]["Domain"]["CDD"]})
    result = get_function(
        "my_id",
        attributes_dict,
        eggnog_annot,
        ipr_test,
        ipr_memberdb_only,
        GeneCaller.PROKKA,
    )
    assert result == ("Pfam Fam function", "InterPro(Pfam)")


def test_get_function_level_pfam_overtake_rejected(setup_data):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    ipr_test = copy.deepcopy(ipr_info)
    del ipr_test["my_id"]["Family"]["NCBIFam"]
    ipr_test["my_id"]["Family"].update({"CDD": ipr_test["my_id"]["Domain"]["CDD"]})
    ipr_test["my_id"]["Family"]["CDD"].update({"level": 1})
    result = get_function(
        "my_id",
        attributes_dict,
        eggnog_annot,
        ipr_test,
        ipr_memberdb_only,
        GeneCaller.PROKKA,
    )
    assert result == ("CDD Domain function", "InterPro(CDD)")


def test_insert_product_source(setup_data):
    _, attributes_dict, _, _ = setup_data
    product_source = "UniFIRE"
    result = insert_product_source(attributes_dict, product_source)
    assert result == {
        "ID": "my_id",
        "product": "hypothetical protein",
        "product_source": "UniFIRE",
        "protein_id": "protein id",
    }


def test_escape_reserved_characters():
    func = "2,3-containing domain, functional protein; integrase"
    result = escape_reserved_characters(func)
    assert result == "2%2C3-containing domain/ functional protein/ integrase"


def test_clean_up_function_domain():
    func = "ATPase, nucleotide binding domain"
    result = clean_up_function(func)
    assert result == "ATPase, nucleotide binding domain-containing protein"


def test_clean_up_function_no_change():
    func = "Ribosomal protein L23/L15e core domain superfamily"
    result = clean_up_function(func)
    assert result == func


def test_clean_up_function_add_protein():
    func = "Ribosomal L23/L15e core domain superfamily"
    result = clean_up_function(func)
    assert result == "Ribosomal L23/L15e core domain superfamily protein"


def test_keep_or_move_to_note_sentence():
    func = "This is a protein"
    result = keep_or_move_to_note(
        func, "eggNOG", {"ID": "123", "locus_tag": "locus_tag"}, GeneCaller.PROKKA
    )
    assert result == (
        "hypothetical protein",
        "Prokka",
        {"ID": "123", "locus_tag": "locus_tag", "note": "eggNOG:This is a protein"},
    )


def test_keep_or_move_to_note_correct():
    func = "ABC-hydrolase interfering protein with extra large domain"
    result = keep_or_move_to_note(
        func, "eggNOG", {"ID": "123", "locus_tag": "locus_tag"}, GeneCaller.PROKKA
    )
    assert result == (
        func,
        "eggNOG",
        {"ID": "123", "locus_tag": "locus_tag"},
    )


def test_move_function_to_note_with_existing_note():
    col9_with_note = ";".join(
        [
            "ID=BU_ATCC8492_03165",
            "locus_tag=BU_ATCC8492_03165",
            "note=UPF0056 inner membrane protein YhgN",
            "product=NAAT family transporter",
            "product_source=NCBIfam",
            "eggNOG=585543.HMPREF0969_01099",
        ]
    )

    col9_dict = dict(
        re.split(r"(?<!\\)=", item) for item in re.split(r"(?<!\\);", col9_with_note)
    )

    expected_result = {
        "ID": "BU_ATCC8492_03165",
        "locus_tag": "BU_ATCC8492_03165",
        "note": "UPF0056 inner membrane protein YhgN, eggNOG:some function",
        "product": "NAAT family transporter",
        "product_source": "NCBIfam",
        "eggNOG": "585543.HMPREF0969_01099",
    }

    result = move_function_to_note("some function", col9_dict)
    assert result == expected_result


def test_move_function_to_note_with_no_note():
    col9_no_note = ";".join(
        [
            "ID=BU_ATCC8492_00043",
            "inference=ab initio prediction:Prodigal:002006",
            "locus_tag=BU_ATCC8492_00043",
            "product=Protein of unknown function DUF2961",
            "product_source=InterPro(Pfam)",
            "eggNOG=411479.BACUNI_03969",
        ]
    )
    col9_dict = dict(
        re.split(r"(?<!\\)=", item) for item in re.split(r"(?<!\\);", col9_no_note)
    )
    expected_result = {
        "ID": "BU_ATCC8492_00043",
        "inference": "ab initio prediction:Prodigal:002006",
        "locus_tag": "BU_ATCC8492_00043",
        "note": "eggNOG:some function",
        "product": "Protein of unknown function DUF2961",
        "product_source": "InterPro(Pfam)",
        "eggNOG": "411479.BACUNI_03969",
    }
    result = move_function_to_note("some function", col9_dict)
    assert result == expected_result
