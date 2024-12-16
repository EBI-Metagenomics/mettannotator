process PSEUDOFINDER {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/pseudofinder:1.1.0'

    input:
    tuple val(meta), path(compliant_gbk)
    tuple path(pseudofinder_db), val(db_version)

    output:
    tuple val(meta), file("${meta.prefix}_pseudos.gff"), emit: pseudofinder_gff
    path "versions.yml" , emit: versions

    script:
    """
    pseudofinder.py annotate \
    -g ${compliant_gbk} \
    -db ${pseudofinder_db} \
    -op ${meta.prefix} \
    -t ${task.cpus} \
    --diamond

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pseudofinder: \$(python -c "import pseudofinder; print(pseudofinder.__version__)")
        Swiss-Prot: $db_version
    END_VERSIONS
    """
}

process PSEUDOFINDER_POSTPROCESSING {

    tag "${meta.prefix}"

    label 'process_nano'

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta), path(annotations_gff), path(compliant_gff, stageAs: "compliant/*"), path(pseudofinder_gff)

    output:
    tuple val(meta), file("${meta.prefix}_processed_pseudogenes.gff"), emit: pseudofinder_processed_gff
    path "versions.yml" , emit: versions

    script:
    """
    adjust_pseudofinder_output.py \
    --pseudofinder-output ${pseudofinder_gff} \
    --standard-gff ${annotations_gff} \
    --compliant-gff ${compliant_gff} \
    -o ${meta.prefix}_processed_pseudogenes.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
