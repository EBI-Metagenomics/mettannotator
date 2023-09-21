process AMRFINDER_PLUS {

    tag "${meta.prefix}"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.11.4--h6e70893_0'

    input:
    tuple val(meta), path(fna), path(faa), path(gff)

    output:
    tuple val(meta), path("${meta.prefix}_amrfinderplus.tsv"), emit: amrfinder_tsv

    script:
    """
    amrfinder --plus \
    -n ${fna} \
    -p ${faa} \
    -g ${gff} \
    -d ${params.amrfinder_plus_db} \
    -a prokka \
    --output ${meta.prefix}_amrfinderplus.tsv \
    --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(amrfinder --version)
    END_VERSIONS
    """
}
