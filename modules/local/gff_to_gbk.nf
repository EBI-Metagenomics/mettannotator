process GFF_TO_GBK {

    tag "${meta.prefix}"

    label 'process_single'

    container 'quay.io/microbiome-informatics/mgnify-pipelines-toolkit:1.4.19'

    input:
    tuple val(meta), path(gff), path(faa), path(contigs)

    output:
    tuple val(meta), path("${meta.prefix}.gbk"), emit: gbk
    path "versions.yml",                          emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    gbk_generator \\
        --contigs ${contigs} \\
        --gff ${gff} \\
        --faa ${faa} \\
        --output_gbk ${meta.prefix}.gbk \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.prefix}.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
