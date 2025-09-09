process ANTISMASH {

    tag "${meta.prefix}"

    container "nf-core/antismash:8.0.1--pyhdfd78af_0"

    input:
    tuple val(meta), path(gbk)
    tuple path(antismash_db), val(db_version)

    output:
    tuple val(meta), path("${meta.prefix}_results/${meta.prefix}.gbk"), emit: gbk
    tuple val(meta), path("${meta.prefix}_antismash.tar.gz")          , emit: results_tarball
    tuple val(meta), path("${meta.prefix}_antismash.gff")             , emit: gff
    path "versions.yml"                                               , emit: versions

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

    antismash_to_gff.py \\
        -r ${meta.prefix}_results/${meta.prefix}.json -a \$(echo \$(antismash --version | sed 's/^antiSMASH //' )) \\
        -o ${meta.prefix}_antismash.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
        antiSMASH database: $db_version
    END_VERSIONS
    """
}
