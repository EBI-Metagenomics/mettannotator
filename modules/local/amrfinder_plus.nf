process AMRFINDER_PLUS {

    tag "${meta.prefix}"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'

    input:
    tuple val(meta), path(fna), path(faa), path(gff)
    tuple path(amrfinder_plus_db), val(db_version)

    output:
    tuple val(meta), path("${meta.prefix}_amrfinderplus.tsv"), emit: amrfinder_tsv
    path "versions.yml"                                      , emit: versions

    script:
    """
    amrfinder --plus \
    -n ${fna} \
    -p ${faa} \
    -g ${gff} \
    -d ${amrfinder_plus_db} \
    -a prokka \
    --output ${meta.prefix}_amrfinderplus.tsv \
    --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amdfinderplus database: $db_version
    END_VERSIONS
    """
}

process AMRFINDER_PLUS_TO_GFF {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta),path(amrfinder_tsv)

    output:
    tuple val(meta), path("${meta.prefix}_amrfinderplus.gff"), emit: amrfinder_gff
    path "versions.yml", emit: versions

    script:
    // TODO: AMRFINDER_PLUS and AMRFINDER_PLUS_TO_GFF should be just one module
    //       that requires AMRFINDER_PLUS container to be modified and python included
    //       this is the reason for the hardcoded version parameter (-v)
    """
    process_amrfinderplus_results.py \\
    -i ${amrfinder_tsv} \\
    -o ${meta.prefix}_amrfinderplus.gff \\
    -v 3.12.8

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
