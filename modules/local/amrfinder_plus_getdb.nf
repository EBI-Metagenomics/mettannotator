
process AMRFINDER_PLUS_GETDB {

    tag "AMR Finder DB 2026-05-15.1"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.2.7--hf69ffd2_0' :
        'biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0' }"

    // amrfinder_index won't work if singularity mounts the workdir
    scratch false

    publishDir "$params.dbs", mode: 'copy'

    output:
    tuple path("amrfinder", type: "dir"), val("2026-05-15.1"), emit: amrfinder_plus_db

    script:
    """
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$PWD/tmp"

    wget -r -nH --cut-dirs=5 \\
    ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/2026-05-15.1/

    amrfinder_index 2026-05-15.1

    mv 2026-05-15.1 amrfinder

    echo "2026-05-15.1" > amrfinder/VERSION.txt
    """
}
