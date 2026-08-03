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
    # pseudofinder --diamond requires a pre-built <fasta>.dmnd alongside the FASTA and does
    # not build it itself; the staged DB ships only uniprot_sprot.fasta, so build the index
    # here in the (writable) work dir using the same diamond that runs the search.
    diamond makedb \
    --in ${pseudofinder_db}/uniprot_sprot.fasta \
    --db uniprot_sprot.fasta \
    --threads ${task.cpus}

    pseudofinder.py annotate \
    -g ${compliant_gbk} \
    -db uniprot_sprot.fasta \
    -op ${meta.prefix} \
    -t ${task.cpus} \
    --diamond

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pseudofinder: 1.1.0
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
