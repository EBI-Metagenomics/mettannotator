#!/bin/env python3

import copy
import os
import re
import sys
import unittest

sys.path.append(os.path.sep.join(os.getcwd().split(os.path.sep)[:-2]))

from bin.add_hypothetical_protein_descriptions import (
    keep_or_move_to_note,
    move_function_to_note,
    clean_up_function,
    insert_product_source,
    get_function,
    escape_reserved_characters,
)


class TestHypotheticalProteinAnnotation(unittest.TestCase):
    def setUp(self):
        self.eggnog_annot = {"my_id": "eggnog function"}
        self.attributes_dict = {
            "ID": "my_id",
            "product": "hypothetical protein",
            "protein_id": "protein id",
        }
        self.ipr_info = {
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
        self.ipr_memberdb_only = {
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

    def test_get_function_ncbifam_priority(self):
        # Tests which function source is best available. In this test, NCBIFam has priority.
        result = get_function(
            "my_id",
            self.attributes_dict,
            self.eggnog_annot,
            self.ipr_info,
            self.ipr_memberdb_only,
        )
        self.assertEqual(result, ("NCBI Fam function", "InterPro(NCBIFam)"))

    def test_get_function_level_priority(self):
        # Tests which function source is the best one available. In this test, Pfam has priority.
        ipr_test = copy.deepcopy(self.ipr_info)
        ipr_test["my_id"]["Family"]["Pfam"].update({"level": 1})
        result = get_function(
            "my_id",
            self.attributes_dict,
            self.eggnog_annot,
            ipr_test,
            self.ipr_memberdb_only,
        )
        self.assertEqual(result, ("Pfam Fam function", "InterPro(Pfam)"))

    def test_get_function_level_eggnog(self):
        # Tests which function source is the best one available. In this test, eggnog is the only result.
        result = get_function(
            "my_id",
            self.attributes_dict,
            self.eggnog_annot,
            {},
            {},
        )
        self.assertEqual(result, ("eggnog function", "eggNOG"))

    def test_get_function_level_pfam_overtake(self):
        # Tests which function source is the best one available. In this test, Pfam overtakes other annotations
        # because it has priority.
        ipr_test = copy.deepcopy(self.ipr_info)
        del ipr_test["my_id"]["Family"]["NCBIFam"]
        ipr_test["my_id"]["Family"].update({"CDD": ipr_test["my_id"]["Domain"]["CDD"]})
        result = get_function(
            "my_id",
            self.attributes_dict,
            self.eggnog_annot,
            ipr_test,
            self.ipr_memberdb_only,
        )
        self.assertEqual(result, ("Pfam Fam function", "InterPro(Pfam)"))

    def test_get_function_level_pfam_overtake_rejected(self):
        # Tests which function source is the best one available. In this test, Pfam should not take over other
        # annotations because it has a higher level.
        ipr_test = copy.deepcopy(self.ipr_info)
        del ipr_test["my_id"]["Family"]["NCBIFam"]
        ipr_test["my_id"]["Family"].update({"CDD": ipr_test["my_id"]["Domain"]["CDD"]})
        ipr_test["my_id"]["Family"]["CDD"].update({"level": 1})
        result = get_function(
            "my_id",
            self.attributes_dict,
            self.eggnog_annot,
            ipr_test,
            self.ipr_memberdb_only,
        )
        self.assertEqual(result, ("CDD Domain function", "InterPro(CDD)"))

    def test_insert_product_source(self):
        # Tests adding the product source field to the 9th column, after the product field
        product_source = "UniFIRE"
        result = insert_product_source(self.attributes_dict, product_source)
        self.assertEqual(
            result,
            {
                "ID": "my_id",
                "product": "hypothetical protein",
                "product_source": "UniFIRE",
                "protein_id": "protein id",
            },
        )

    def test_escape_reserved_characters(self):
        # Tests replacing special characters
        func = "2,3-containing domain, functional protein; integrase"
        result = escape_reserved_characters(func)
        self.assertEqual(
            result, "2%2C3-containing domain/ functional protein\; integrase"
        )

    def test_clean_up_function_domain(self):
        # Description ends with "domain", need to add "-containing protein"
        func = "ATPase, nucleotide binding domain"
        result = clean_up_function(func)
        self.assertEqual(result, "ATPase, nucleotide binding domain-containing protein")

    def test_clean_up_function_no_change(self):
        # Description contains the word "protein", no need to change
        func = "Ribosomal protein L23/L15e core domain superfamily"
        result = clean_up_function(func)
        self.assertEqual(result, func)

    def test_clean_up_function_add_protein(self):
        # Description ends with "family" and doesn't contain the word protein, "protein" should be added
        func = "Ribosomal L23/L15e core domain superfamily"
        result = clean_up_function(func)
        self.assertEqual(result, "Ribosomal L23/L15e core domain superfamily protein")

    def test_keep_or_move_to_note_sentence(self):
        # Function is a sentence, should be moved to the note field
        func = "This is a protein"
        result = keep_or_move_to_note(
            func, "eggNOG", {"ID": "123", "locus_tag": "locus_tag"}
        )
        self.assertEqual(
            result,
            (
                "hypothetical protein",
                "Prokka",
                {
                    "ID": "123",
                    "locus_tag": "locus_tag",
                    "note": "eggNOG:This is a protein",
                },
            ),
        )

    def test_keep_or_move_to_note_correct(self):
        # Function is in the correct format, keep it as is.
        func = "ABC-hydrolase interfering protein with extra large domain"
        result = keep_or_move_to_note(
            func, "eggNOG", {"ID": "123", "locus_tag": "locus_tag"}
        )
        self.assertEqual(
            result,
            (
                func,
                "eggNOG",
                {"ID": "123", "locus_tag": "locus_tag"},
            ),
        )

    def test_move_function_to_note_with_existing_note(self):
        # Tests appending function description to the existing text in the note field
        col9_with_note = "ID=BU_ATCC8492_03165;locus_tag=BU_ATCC8492_03165;note=UPF0056 inner membrane protein YhgN;product=NAAT family transporter;product_source=NCBIfam;eggNOG=585543.HMPREF0969_01099"
        col9_dict = dict(
            re.split(r"(?<!\\)=", item)
            for item in re.split(r"(?<!\\);", col9_with_note)
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
        self.assertEqual(result, expected_result)

    def test_move_function_to_note_with_no_note(self):
        # Tests creating a note field where it doesn't already exist and saving the function into it
        col9_no_note = "ID=BU_ATCC8492_00043;inference=ab initio prediction:Prodigal:002006;locus_tag=BU_ATCC8492_00043;product=Protein of unknown function DUF2961;product_source=InterPro(Pfam);eggNOG=411479.BACUNI_03969"
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
        self.assertEqual(result, expected_result)


if __name__ == "__main__":
    unittest.main()
