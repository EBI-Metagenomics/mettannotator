import pytest

from bin.normalize_ensembl_gff import (
    col9_to_dict,
    decode_description,
    normalise,
    parse_col9,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

GENE = (
    "contig_1\tEnsembl\tgene\t100\t300\t.\t+\t.\t"
    "ID=gene:ABAYE0001;gene_id=ABAYE0001;Name=recA;description=DNA%20recombinase%20A\n"
)
MRNA = (
    "contig_1\tEnsembl\tmRNA\t100\t300\t.\t+\t.\t"
    "ID=transcript:T001;Parent=gene:ABAYE0001;transcript_id=T001\n"
)
CDS_WITH_PROTEIN_ID = (
    "contig_1\tEnsembl\tCDS\t100\t300\t.\t+\t0\t"
    "ID=CDS:PROT001;Parent=transcript:T001;protein_id=PROT001\n"
)
CDS_WITHOUT_PROTEIN_ID = (
    "contig_1\tEnsembl\tCDS\t100\t300\t.\t+\t0\t"
    "ID=CDS:PROT001;Parent=transcript:T001\n"
)


def run(tmp_path, *lines):
    in_f = tmp_path / "in.gff"
    out_f = tmp_path / "out.gff"
    in_f.write_text("##gff-version 3\n" + "".join(lines))
    normalise(str(in_f), str(out_f))
    return out_f.read_text().splitlines()


def cds(lines):
    return [line for line in lines if "\tCDS\t" in line]


def genes(lines):
    return [line for line in lines if "\tgene\t" in line]


# ---------------------------------------------------------------------------
# parse_col9
# ---------------------------------------------------------------------------


def test_parse_col9_basic():
    assert parse_col9("ID=gene:G1;Name=recA;gene_id=G1") == [
        ("ID", "gene:G1"),
        ("Name", "recA"),
        ("gene_id", "G1"),
    ]


def test_parse_col9_empty_segments_ignored():
    assert parse_col9("ID=foo;;Name=bar") == [("ID", "foo"), ("Name", "bar")]


def test_parse_col9_attribute_without_value():
    assert parse_col9("ID=foo;flag") == [("ID", "foo"), ("flag", "")]


def test_parse_col9_value_containing_equals():
    # Only the first = splits key from value
    assert parse_col9("ID=foo;description=a=b") == [
        ("ID", "foo"),
        ("description", "a=b"),
    ]


def test_parse_col9_single_attribute():
    assert parse_col9("ID=only") == [("ID", "only")]


# ---------------------------------------------------------------------------
# col9_to_dict
# ---------------------------------------------------------------------------


def test_col9_to_dict_basic():
    result = col9_to_dict("ID=G1;Name=recA")
    assert result == {"ID": "G1", "Name": "recA"}


def test_col9_to_dict_last_wins_on_duplicate_keys():
    # dict comprehension keeps last value for duplicate keys
    result = col9_to_dict("ID=first;ID=second")
    assert result["ID"] == "second"


# ---------------------------------------------------------------------------
# decode_description
# ---------------------------------------------------------------------------


def test_decode_description_url_encoded_spaces():
    assert decode_description("DNA%20repair%20protein") == "DNA repair protein"


def test_decode_description_semicolons_replaced_with_comma():
    assert decode_description("func1;func2") == "func1,func2"


def test_decode_description_percent_encoded_semicolon():
    assert decode_description("func1%3Bfunc2") == "func1,func2"


def test_decode_description_strips_whitespace():
    assert decode_description("  protein  ") == "protein"


def test_decode_description_empty_string():
    assert decode_description("") == ""


def test_decode_description_no_encoding():
    assert decode_description("kinase") == "kinase"


# ---------------------------------------------------------------------------
# normalise — header handling
# ---------------------------------------------------------------------------


def test_gff_version_directive_kept(tmp_path):
    lines = run(tmp_path)
    assert "##gff-version 3" in lines


def test_sequence_region_directive_kept(tmp_path):
    in_f = tmp_path / "in.gff"
    out_f = tmp_path / "out.gff"
    in_f.write_text("##gff-version 3\n##sequence-region contig_1 1 5000\n")
    normalise(str(in_f), str(out_f))
    assert "##sequence-region contig_1 1 5000" in out_f.read_text().splitlines()


def test_ensembl_pragma_dropped(tmp_path):
    lines = run(tmp_path, "#!genome-build GCA_000007405.1\n")
    assert not any(line.startswith("#!") for line in lines)


def test_separator_lines_dropped(tmp_path):
    lines = run(tmp_path, "###\n###\n")
    assert "###" not in lines


# ---------------------------------------------------------------------------
# normalise — feature filtering
# ---------------------------------------------------------------------------


def test_mrna_feature_dropped(tmp_path):
    lines = run(tmp_path, MRNA)
    assert not any("\tmRNA\t" in line for line in lines)


@pytest.mark.parametrize(
    "feature",
    ["exon", "pseudogene", "ncRNA_gene", "rRNA", "tRNA", "biological_region"],
)
def test_skip_features_dropped(tmp_path, feature):
    feat_line = f"contig_1\tEnsembl\t{feature}\t1\t100\t.\t+\t.\tID=x\n"
    lines = run(tmp_path, feat_line)
    assert not any(f"\t{feature}\t" in line for line in lines)


# ---------------------------------------------------------------------------
# normalise — gene features
# ---------------------------------------------------------------------------


def test_gene_id_gets_gene_suffix(tmp_path):
    gl = genes(run(tmp_path, GENE))
    assert len(gl) == 1
    assert "ID=ABAYE0001_gene" in gl[0]


def test_gene_locus_tag_set(tmp_path):
    gl = genes(run(tmp_path, GENE))
    assert "locus_tag=ABAYE0001" in gl[0]


def test_gene_name_kept_when_present(tmp_path):
    gl = genes(run(tmp_path, GENE))
    assert "Name=recA" in gl[0]


def test_gene_name_omitted_when_absent(tmp_path):
    gene_no_name = (
        "contig_1\tEnsembl\tgene\t1\t100\t.\t+\t.\t"
        "ID=gene:G2;gene_id=G2;description=some%20protein\n"
    )
    gl = genes(run(tmp_path, gene_no_name))
    assert "Name=" not in gl[0]


# ---------------------------------------------------------------------------
# normalise — CDS features
# ---------------------------------------------------------------------------


def test_cds_protein_id_used_as_id(tmp_path):
    cl = cds(run(tmp_path, GENE, MRNA, CDS_WITH_PROTEIN_ID))
    assert len(cl) == 1
    assert "ID=PROT001" in cl[0]
    assert "ID=CDS:" not in cl[0]


def test_cds_prefix_stripped_when_no_protein_id(tmp_path):
    cl = cds(run(tmp_path, GENE, MRNA, CDS_WITHOUT_PROTEIN_ID))
    assert "ID=PROT001" in cl[0]
    assert "ID=CDS:" not in cl[0]


def test_cds_product_from_parent_gene(tmp_path):
    cl = cds(run(tmp_path, GENE, MRNA, CDS_WITH_PROTEIN_ID))
    assert "product=DNA recombinase A" in cl[0]


def test_cds_product_hypothetical_when_no_gene_info(tmp_path):
    cl = cds(run(tmp_path, CDS_WITH_PROTEIN_ID))
    assert "product=hypothetical protein" in cl[0]


def test_cds_has_no_parent_attribute(tmp_path):
    cl = cds(run(tmp_path, GENE, MRNA, CDS_WITH_PROTEIN_ID))
    assert "Parent=" not in cl[0]


def test_cds_locus_tag_from_parent_gene(tmp_path):
    cl = cds(run(tmp_path, GENE, MRNA, CDS_WITH_PROTEIN_ID))
    assert "locus_tag=ABAYE0001" in cl[0]


def test_cds_description_url_decoded(tmp_path):
    gene = (
        "contig_1\tEnsembl\tgene\t1\t100\t.\t+\t.\t"
        "ID=gene:G3;gene_id=G3;description=DNA%20repair%20protein\n"
    )
    mrna = (
        "contig_1\tEnsembl\tmRNA\t1\t100\t.\t+\t.\t"
        "ID=transcript:T003;Parent=gene:G3;transcript_id=T003\n"
    )
    cds_line = (
        "contig_1\tEnsembl\tCDS\t1\t100\t.\t+\t0\t"
        "ID=CDS:P003;Parent=transcript:T003;protein_id=P003\n"
    )
    cl = cds(run(tmp_path, gene, mrna, cds_line))
    assert "product=DNA repair protein" in cl[0]


def test_cds_description_semicolons_replaced(tmp_path):
    gene = (
        "contig_1\tEnsembl\tgene\t1\t100\t.\t+\t.\t"
        "ID=gene:G4;gene_id=G4;description=func1%3Bfunc2\n"
    )
    mrna = (
        "contig_1\tEnsembl\tmRNA\t1\t100\t.\t+\t.\t"
        "ID=transcript:T004;Parent=gene:G4;transcript_id=T004\n"
    )
    cds_line = (
        "contig_1\tEnsembl\tCDS\t1\t100\t.\t+\t0\t"
        "ID=CDS:P004;Parent=transcript:T004;protein_id=P004\n"
    )
    cl = cds(run(tmp_path, gene, mrna, cds_line))
    assert "product=func1,func2" in cl[0]
