process AMRFINDER_PLUS_GET_SPECIES_NAME {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta), path(fasta)
    tuple path(amrfinder_plus_db), val(db_version)

    output:
    tuple val(meta), env(detected_organism), emit: detected_organism

    script:
    // Script checks available organisms in AMRFinderPlus and uses user-provided taxid to determine the organism
    // name to use with the -O option in AMRFinderPlus. If there is no matching organism (taxonomy is not species-level
    // or species is not in the AMRFinderPlus list, script returns an empty string)
    """
    detected_organism=\$(get_species_name_for_amrfinderplus.py -t ${meta.taxid} -d ${amrfinder_plus_db})
    """

    stub:
    """
    detected_organism=""
    """

}

process AMRFINDER_PLUS {

    tag "${meta.prefix}"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.23--hf69ffd2_0'

    input:
    tuple val(meta), path(fna), path(faa), path(gff), val(detected_organism)
    tuple path(amrfinder_plus_db), val(db_version)

    output:
    tuple val(meta), path("${meta.prefix}_amrfinderplus.tsv"), emit: amrfinder_tsv
    tuple val(meta), path("organism.txt")                    , emit: amrfinderplus_organism_file
    path "versions.yml"                                      , emit: versions

    script:
    """
    # this is needed as some environments are very picky with the TMP dir
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$PWD/tmp"
    """

    if ( detected_organism == "" )
        organism_flag = ""
        """
        echo "--organism option not used - taxon is not present in the AMRFinderPlus list" > organism.txt
        """
    else
        organism_flag = "-O ${detected_organism}"
        """
        echo \$detected_organism > organism.txt
        """
    """
    amrfinder --plus \
    -n ${fna} \
    -p ${faa} \
    -g ${gff} \
    -d ${amrfinder_plus_db} \
    -a prokka \
    --output ${meta.prefix}_amrfinderplus.tsv \
    --threads ${task.cpus} ${organism_flag}

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
    -v 4.0.23

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
