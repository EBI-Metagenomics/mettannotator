process ANTISMASH {

    tag "${meta.prefix}"

    label "error_retry"

    container "nf-core/antismash:8.0.1--pyhdfd78af_0"

    input:
    tuple val(meta), path(gbk)
    tuple path(antismash_db), val(db_version)

    output:
    tuple val(meta), path("${meta.prefix}_results/${meta.prefix}.gbk") , emit: gbk
    tuple val(meta), path("${meta.prefix}_antismash.tar.gz")           , emit: results_tarball
    tuple val(meta), path("${meta.prefix}_results/${meta.prefix}.json"), emit: json
    path "versions.yml"                                                , emit: versions

    script:
    """
    antismash \\
    -t bacteria \\
    -c ${task.cpus} \\
    --databases ${antismash_db} \\
    --output-basename ${meta.prefix} \\
    --genefinding-tool none \\
    --allow-long-headers \\
    --output-dir ${meta.prefix}_results \\
    ${gbk}

    tar -czf ${meta.prefix}_antismash.tar.gz ${meta.prefix}_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
        antiSMASH database: $db_version
    END_VERSIONS
    """
}

process ANTISMASH_TO_GFF {

    tag "${meta.prefix}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.2.9--pyhdfd78af_0' :
        'biocontainers/mgnify-pipelines-toolkit:1.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(antismash_json)

    output:
    tuple val(meta), path("${meta.prefix}_antismash.gff"), emit: gff
    path "versions.yml", emit: versions

    script:
    """
    antismash_gff_builder \\
    -i ${antismash_json} \\
    -o ${meta.prefix}_antismash.gff \\
    --cds_tag locus_tag

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}

process ANTISMASH_SUMMARY {

    tag "${meta.prefix}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.2.9--pyhdfd78af_0' :
        'biocontainers/mgnify-pipelines-toolkit:1.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(antismash_gff)

    output:
    tuple val(meta), path("*_summary.tsv")       , emit: antismash_summary
    path "versions.yml"                          , emit: versions

    script:
    """
    summarise_antismash_bgcs \\
        --antismash-gff ${antismash_gff} \\
        --output ${antismash_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

}
