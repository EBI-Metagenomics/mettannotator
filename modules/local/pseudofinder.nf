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
    -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pseudofinder: \$(python -c "import pseudofinder; print(pseudofinder.__version__)")
        Swiss-Prot: $db_version
    END_VERSIONS
    """
}
