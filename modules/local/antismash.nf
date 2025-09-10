process ANTISMASH {

    tag "${meta.prefix}"

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

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta), path(antismash_json)

    output:
    tuple val(meta), path("${meta.prefix}_antismash.gff"), emit: gff
    path "versions.yml", emit: versions

    script:
    """
    antismash_gff_builder.py \\
    -i ${antismash_json} \\
    -o ${meta.prefix}_antismash.gff \\
    --cds_tag locus_tag

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
