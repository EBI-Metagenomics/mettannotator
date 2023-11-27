process ANTISMASH {

    tag "${meta.prefix}"

    label "process_low"

    container 'quay.io/microbiome-informatics/antismash:7.1.0.1'

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${meta.prefix}_results/${meta.prefix}.gbk"), emit: gbk
    tuple val(meta), path("${meta.prefix}_antismash.tar.gz")          , emit: results_tarball
    path "versions.yml"                                               , emit: versions

    script:
    """
    antismash \
    -t bacteria \
    -c ${task.cpus} \
    --databases ${params.antismash_db} \
    --output-basename ${meta.prefix} \
    --output-dir ${meta.prefix}_results \
    ${gbk}

    tar -czf ${meta.prefix}_antismash.tar.gz ${meta.prefix}_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
    END_VERSIONS
    """
}
