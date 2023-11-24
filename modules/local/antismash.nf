process ANTISMASH {

    tag "${meta.prefix}"

    label "process_low"

    container 'quay.io/microbiome-informatics/antismash:7.1.0.1'

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${meta.prefix}_results/*"), emit: bgc_results
    path "versions.yml", emit: versions

    script:
    """
    antismash \
    -t bacteria \
    -c ${task.cpus} \ 
    --databases ${params.antismash_db} \
    --output-basename ${meta.prefix} \
    --output-dir ${meta.prefix}_results \
    ${gbk}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
    END_VERSIONS
    """
}
