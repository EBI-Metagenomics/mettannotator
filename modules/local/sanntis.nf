/*
 * SMBGC Annotation using Neural Networks Trained on Interpro Signatures
*/
process SANNTIS {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/sanntis:0.9.3.2'

    input:
    tuple val(meta), path(interproscan_tsv), path(prokka_gbk)

    output:
    tuple val(meta), path("*_sanntis.gff"), emit: sanntis_gff
    path "versions.yml"                   , emit: versions

    script:
    if (interproscan_tsv.extension == "gz") {
        """
        gunzip -c ${interproscan_tsv} > interproscan.tsv
        sanntis \
        --ip-file interproscan.tsv \
        --outfile ${meta.prefix}_sanntis.gff \
        ${prokka_gbk}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
        END_VERSIONS
        """
    } else {
        """
        sanntis \
        --ip-file ${interproscan_tsv} \
        --outfile ${meta.prefix}_sanntis.gff \
        ${prokka_gbk}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
        END_VERSIONS
        """
    }
}
