/*
 * SMBGC Annotation using Neural Networks Trained on Interpro Signatures
*/
process SANNTIS {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/sanntis:0.9.3.4'

    input:
    tuple val(meta), path(interproscan_tsv), path(prokka_gbk)

    output:
    tuple val(meta), path("*_sanntis.gff"), emit: sanntis_gff
    path "versions.yml"                   , emit: versions

    script:
    def genbank_extension = prokka_gbk.extension
    if (interproscan_tsv.extension == "gz") {
        """
        gunzip -c ${interproscan_tsv} > interproscan.tsv
        // Workaround implemented for SanntiS due to discrepancies in Bakta's format, resulting in empty output files.
        grep -v "/protein_id=" ${prokka_gbk} > ${meta.prefix}_prepped.${genbank_extension}
        sanntis \
        --ip-file interproscan.tsv \
        --outfile ${meta.prefix}_sanntis.gff \
        --cpu ${task.cpus} \
        ${meta.prefix}_prepped.${genbank_extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
        END_VERSIONS
        """
    } else {
        """
        grep -v "/protein_id=" ${prokka_gbk} > ${meta.prefix}_prepped.${genbank_extension}
        sanntis \
        --ip-file ${interproscan_tsv} \
        --outfile ${meta.prefix}_sanntis.gff \
        ${meta.prefix}_prepped.${genbank_extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
        END_VERSIONS
        """
    }
}
