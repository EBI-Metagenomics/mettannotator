process AMRFINDER_PLUS {

    tag "${meta.prefix}"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.11.4--h6e70893_0'

    input:
    tuple val(meta), path(fna), path(faa), path(gff)

    output:
    tuple val(meta), path("${meta.prefix}_amrfinderplus.tsv"), emit: amrfinder_tsv
    tuple val(meta), path("${meta.prefix}_amrfinderplus.gff"), emit: amrfinder_gff
    path "versions.yml", emit: versions

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

    process_amrfinderplus_results.py -i ${meta.prefix}_amrfinderplus.tsv -o ${meta.prefix}_amrfinderplus.gff -v 3.11.4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}
