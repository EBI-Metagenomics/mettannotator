process RENAME_CONTIGS {
    tag "$meta.prefix"
    label 'process_low'

    conda "bioconda::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(input_files)
    path mapping_file, stageAs: "mapping.tsv"

    output:
    tuple val(meta), path("renamed/*.{fa,fasta,fna}"), emit: modified_fasta, optional: true
    tuple val(meta), path("renamed/*.{gbk}")         , emit: modified_gbk, optional: true
    tuple val(meta), path("renamed/*.{gff}")         , emit: modified_gff, optional: true
    tuple val(meta), path("*.map.tsv")               , emit: fasta_ids_mapping, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mapping_arg = mapping_file ? "--map ${mapping_file}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p renamed

    rename_contigs.py \\
        ${input_files} \\
        --outdir renamed \\
        ${mapping_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p renamed

    # Create stub output files based on input file extensions
    for file in ${input_files}; do
        touch "renamed/\${file}"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
