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
    UnirefHit,
    add_uniref_attributes,
    clean_up_function,
    escape_reserved_characters,
    get_function,
    insert_product_source,
    is_informative_uniref_name,
    keep_or_move_to_note,
    load_uniref,
    main,
    move_function_to_note,
    parse_uniref_title,
)

DIAMOND_HEADER = (
    "qseqid\tstitle\tqlen\tslen\tlength\tqcovhsp\tpident\tevalue\tbitscore\tscovhsp"
)


def diamond_line(
    acc, cluster, name, pident, evalue="1e-50", bitscore="300", qcov="100", scov="95.2"
):
    rep_id = cluster.split("_", 1)[1]
    stitle = f"{cluster} {name} n=3 Tax=Bacteroides TaxID=816 RepID={rep_id}_BACXX"
    return "\t".join(
        [acc, stitle, "200", "210", "200", qcov, str(pident), evalue, bitscore, scov]
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
        "prokka",
        {"ID": "123", "locus_tag": "locus_tag", "Note": "eggNOG:This is a protein"},
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
            "Note=UPF0056 inner membrane protein YhgN",
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
        "Note": "UPF0056 inner membrane protein YhgN, eggNOG:some function",
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
        "Note": "eggNOG:some function",
        "product": "Protein of unknown function DUF2961",
        "product_source": "InterPro(Pfam)",
        "eggNOG": "411479.BACUNI_03969",
    }
    result = move_function_to_note("some function", col9_dict)
    assert result == expected_result


# ---------------------------------------------------------------------------
# UniRef
# ---------------------------------------------------------------------------


def test_parse_uniref_title():
    stitle = (
        "UniRef90_A0A849BTK1 NAD(+)--protein-arginine ADP-ribosyltransferase n=1 "
        "Tax=Nocardia uniformis TaxID=53432 RepID=A0A849BTK1_9NOCA"
    )
    assert parse_uniref_title(stitle) == (
        "UniRef90_A0A849BTK1",
        "NAD(+)--protein-arginine ADP-ribosyltransferase",
    )


def test_parse_uniref_title_strips_fragment():
    stitle = (
        "UniRef50_A0A949BZA8 Adenylate kinase (Fragment) n=2 "
        "Tax=Omnitrophica bacterium TaxID=2035772 RepID=A0A949BZA8_UNCOM"
    )
    assert parse_uniref_title(stitle) == ("UniRef50_A0A949BZA8", "Adenylate kinase")


def test_parse_uniref_title_not_uniref():
    assert parse_uniref_title("sp|P0A7B8|KAD_ECOLI Adenylate kinase") is None


@pytest.mark.parametrize(
    "name",
    [
        "Uncharacterized protein",
        "uncharacterized protein",
        "Putative uncharacterized protein",
        "Uncharacterised protein",
        "Uncharacterized protein YhgN",
        "Uncharacterized conserved protein",
        "Uncharacterized protein conserved in bacteria",
        "Uncharacterized oxidoreductase YdgJ",
        "Hypothetical protein",
        "hypothetical protein",
        "Conserved hypothetical protein",
        "Predicted hypothetical protein",
        "Protein of unknown function",
        "Protein of unknown function DUF1232",
        "Domain of unknown function",
        "Predicted protein",
        "Unknown protein",
        "DUF4974 domain-containing protein",
        "DUF1016 family protein",
        "",
    ],
)
def test_uninformative_uniref_names(name):
    assert not is_informative_uniref_name(name)


@pytest.mark.parametrize(
    "name",
    [
        "Adenylate kinase",
        "Ig-like domain-containing protein",
        "SusC/RagA family TonB-linked outer membrane protein",
        "UPF0056 inner membrane protein YhgN",
        "Putative ABC transporter permease",
    ],
)
def test_informative_uniref_names(name):
    assert is_informative_uniref_name(name)


def test_load_uniref_filters_and_picks_best_informative_hit(tmp_path):
    tsv = tmp_path / "hits.tsv"
    tsv.write_text(
        "\n".join(
            [
                DIAMOND_HEADER,
                # best score is uninformative, so the next named hit is used
                diamond_line(
                    "p1", "UniRef90_AAA", "Uncharacterized protein", 99, bitscore="500"
                ),
                diamond_line(
                    "p1",
                    "UniRef90_BBB",
                    "Adenylate kinase",
                    95,
                    bitscore="450",
                    qcov="98.5",
                    scov="97.1",
                ),
                diamond_line(
                    "p1", "UniRef90_CCC", "Guanylate kinase", 94, bitscore="400"
                ),
                # below the identity cutoff
                diamond_line("p2", "UniRef90_DDD", "Thioredoxin", 49.9),
                # e-value too high
                diamond_line("p3", "UniRef90_EEE", "Thioredoxin", 99, evalue="1e-5"),
                # only uninformative hits
                diamond_line("p4", "UniRef90_FFF", "Hypothetical protein", 99),
                # higher score later in the file wins
                diamond_line("p5", "UniRef90_GGG", "Ferredoxin", 92, bitscore="100"),
                diamond_line("p5", "UniRef90_HHH", "Rubredoxin", 98, bitscore="200"),
            ]
        )
        + "\n"
    )
    assert load_uniref(str(tsv)) == {
        "p1": UnirefHit("Adenylate kinase", "UniRef90_BBB", "95", "98.5", "97.1"),
        "p5": UnirefHit("Rubredoxin", "UniRef90_HHH", "98", "100", "95.2"),
    }


def test_load_uniref_identity_cutoff_is_50_and_inclusive(tmp_path):
    tsv = tmp_path / "hits.tsv"
    tsv.write_text(
        "\n".join(
            [
                diamond_line("p1", "UniRef90_AAA", "Thioredoxin", "50.0"),
                diamond_line("p2", "UniRef90_BBB", "Thioredoxin", "49.9"),
            ]
        )
        + "\n"
    )
    assert load_uniref(str(tsv)) == {
        "p1": UnirefHit("Thioredoxin", "UniRef90_AAA", "50.0", "100", "95.2")
    }


def test_load_uniref_empty_or_missing(tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    assert load_uniref(str(empty)) == {}
    assert load_uniref(None) == {}


@pytest.fixture
def uniref_hits():
    uniref90 = {
        "my_id": UnirefHit("UniRef90 function", "UniRef90_AAA", "95", "100", "99")
    }
    uniref50 = {
        "my_id": UnirefHit("UniRef50 function", "UniRef50_BBB", "60", "90", "85")
    }
    return uniref90, uniref50


def call_get_function(attributes_dict, eggnog, ipr, memberdb, uniref90, uniref50):
    return get_function(
        "my_id",
        attributes_dict,
        eggnog,
        ipr,
        memberdb,
        GeneCaller.PROKKA,
        uniref90,
        uniref50,
    )


def test_get_function_unifire_beats_uniref(setup_data, uniref_hits):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    attributes = dict(attributes_dict, uf_prot_rec_fullname="UniFIRE function")
    result = call_get_function(
        attributes, eggnog_annot, ipr_info, ipr_memberdb_only, *uniref_hits
    )
    assert result == ("UniFIRE function", "UniFIRE")


def test_get_function_interpro_family_beats_uniref90(setup_data, uniref_hits):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    result = call_get_function(
        attributes_dict, eggnog_annot, ipr_info, ipr_memberdb_only, *uniref_hits
    )
    assert result == ("NCBI Fam function", "InterPro(NCBIFam)")


def test_get_function_uniref90_beats_interpro_domain(setup_data, uniref_hits):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    ipr_test = copy.deepcopy(ipr_info)
    del ipr_test["my_id"]["Family"]
    result = call_get_function(
        attributes_dict, eggnog_annot, ipr_test, ipr_memberdb_only, *uniref_hits
    )
    assert result == ("UniRef90 function", "UniRef90")


def test_get_function_interpro_domain_beats_uniref50(setup_data, uniref_hits):
    eggnog_annot, attributes_dict, ipr_info, ipr_memberdb_only = setup_data
    ipr_test = copy.deepcopy(ipr_info)
    del ipr_test["my_id"]["Family"]
    _, uniref50 = uniref_hits
    result = call_get_function(
        attributes_dict, eggnog_annot, ipr_test, ipr_memberdb_only, {}, uniref50
    )
    assert result == ("CDD Domain function", "InterPro(CDD)")


def test_get_function_memberdb_beats_uniref50(setup_data, uniref_hits):
    eggnog_annot, attributes_dict, _, ipr_memberdb_only = setup_data
    _, uniref50 = uniref_hits
    result = call_get_function(
        attributes_dict, eggnog_annot, {}, ipr_memberdb_only, {}, uniref50
    )
    assert result == ("NCBI Fam function", "InterPro(NCBIFam)")


def test_get_function_uniref50_beats_superfamily_and_eggnog(setup_data, uniref_hits):
    eggnog_annot, attributes_dict, _, _ = setup_data
    _, uniref50 = uniref_hits
    ipr_superfamily = {
        "my_id": {
            "Homologous_superfamily": {
                "SUPERFAMILY": {
                    "match": 0.9,
                    "ipr_desc": "Superfamily function",
                    "sig_desc": "sig desc",
                    "level": 0,
                }
            }
        }
    }
    result = call_get_function(
        attributes_dict, eggnog_annot, ipr_superfamily, {}, {}, uniref50
    )
    assert result == ("UniRef50 function", "UniRef50")


def test_get_function_fast_mode_uniref_order(setup_data, uniref_hits):
    # --fast: no InterProScan or UniFIRE, so UniRef90 > UniRef50 > eggNOG
    eggnog_annot, attributes_dict, _, _ = setup_data
    uniref90, uniref50 = uniref_hits
    assert call_get_function(
        attributes_dict, eggnog_annot, {}, {}, uniref90, uniref50
    ) == ("UniRef90 function", "UniRef90")
    assert call_get_function(attributes_dict, eggnog_annot, {}, {}, {}, uniref50) == (
        "UniRef50 function",
        "UniRef50",
    )
    assert call_get_function(attributes_dict, eggnog_annot, {}, {}, {}, {}) == (
        "eggnog function",
        "eggNOG",
    )


def test_main_writes_uniref_product_and_hit_attributes(tmp_path):
    prodigal = "inference=ab initio prediction:Prodigal:002006"
    gff = tmp_path / "in.gff"
    gff.write_text(
        "\n".join(
            [
                "##gff-version 3",
                f"contig_1\tProdigal:002006\tCDS\t1\t300\t.\t+\t0\t"
                f"ID=LOC_00001;{prodigal};locus_tag=LOC_00001;product=hypothetical protein",
                f"contig_1\tProdigal:002006\tCDS\t400\t700\t.\t+\t0\t"
                f"ID=LOC_00002;{prodigal};locus_tag=LOC_00002;product=hypothetical protein",
                f"contig_1\tProdigal:002006\tCDS\t800\t1100\t.\t+\t0\t"
                f"ID=LOC_00003;{prodigal};locus_tag=LOC_00003;product=hypothetical protein",
                "contig_1\tProdigal:002006\tCDS\t1200\t1500\t.\t+\t0\t"
                "ID=LOC_00004;locus_tag=LOC_00004;product=Adenylate kinase",
            ]
        )
        + "\n"
    )
    eggnog = tmp_path / "eggnog.tsv"
    eggnog.write_text("#query\tseed_ortholog\tevalue\n")
    uniref90 = tmp_path / "uniref90.tsv"
    uniref90.write_text(
        "\n".join(
            [
                DIAMOND_HEADER,
                diamond_line(
                    "LOC_00001",
                    "UniRef90_A0A1",
                    "Thioredoxin; reductase",
                    97.3,
                    qcov="99.1",
                    scov="96.4",
                ),
                diamond_line(
                    "LOC_00003", "UniRef90_A0A3", "Uncharacterized protein", 99
                ),
                diamond_line("LOC_00004", "UniRef90_A0A4", "Something else", 97),
            ]
        )
        + "\n"
    )
    uniref50 = tmp_path / "uniref50.tsv"
    uniref50.write_text(
        "\n".join(
            [
                DIAMOND_HEADER,
                diamond_line("LOC_00001", "UniRef50_B0B1", "Thioredoxin-like", 60),
                diamond_line(
                    "LOC_00002",
                    "UniRef50_B0B2",
                    "Ferredoxin",
                    61.5,
                    qcov="88",
                    scov="90.3",
                ),
                diamond_line("LOC_00003", "UniRef50_B0B3", "Hypothetical protein", 70),
            ]
        )
        + "\n"
    )
    out = tmp_path / "out.gff"

    main(
        None,
        None,
        None,
        str(eggnog),
        str(gff),
        str(out),
        GeneCaller.PROKKA,
        str(uniref90),
        str(uniref50),
    )

    lines = [line for line in out.read_text().splitlines() if not line.startswith("#")]
    col9 = [line.split("\t")[8] for line in lines]
    assert col9[0] == (
        f"ID=LOC_00001;{prodigal};locus_tag=LOC_00001;"
        "product=Thioredoxin/ reductase;product_source=UniRef90;Dbxref=UniRef:UniRef90_A0A1;"
        "uniref_identity=97.3;uniref_qcov=99.1;uniref_scov=96.4"
    )
    assert col9[1] == (
        f"ID=LOC_00002;{prodigal};locus_tag=LOC_00002;"
        "product=Ferredoxin;product_source=UniRef50;Dbxref=UniRef:UniRef50_B0B2;"
        "uniref_identity=61.5;uniref_qcov=88;uniref_scov=90.3"
    )
    # only uncharacterised/hypothetical names were found, so nothing is transferred
    assert col9[2] == (
        f"ID=LOC_00003;{prodigal};locus_tag=LOC_00003;"
        "product=hypothetical protein;product_source=prokka"
    )
    # proteins already named by the gene caller are untouched
    assert col9[3] == (
        "ID=LOC_00004;locus_tag=LOC_00004;product=Adenylate kinase;product_source=prokka"
    )


def test_add_uniref_attributes_appends_to_existing_dbxref():
    # e.g. Bakta output, which already has cross-references
    attributes = {
        "ID": "p1",
        "Dbxref": "SO:0001217,UniParc:UPI0001",
        "product": "Adenylate kinase",
        "product_source": "UniRef90",
        "locus_tag": "p1",
    }
    hit = UnirefHit("Adenylate kinase", "UniRef90_AAA", "95", "100", "99")
    assert add_uniref_attributes(attributes, hit) == {
        "ID": "p1",
        "Dbxref": "SO:0001217,UniParc:UPI0001,UniRef:UniRef90_AAA",
        "product": "Adenylate kinase",
        "product_source": "UniRef90",
        "uniref_identity": "95",
        "uniref_qcov": "100",
        "uniref_scov": "99",
        "locus_tag": "p1",
    }


def test_add_uniref_attributes_does_not_duplicate_dbxref():
    attributes = {
        "ID": "p1",
        "Dbxref": "UniRef:UniRef50_BBB,SO:0001217",
        "product": "Ferredoxin",
        "product_source": "UniRef50",
    }
    hit = UnirefHit("Ferredoxin", "UniRef50_BBB", "60", "90", "85")
    result = add_uniref_attributes(attributes, hit)
    assert result["Dbxref"] == "UniRef:UniRef50_BBB,SO:0001217"
    assert result["uniref_scov"] == "85"
