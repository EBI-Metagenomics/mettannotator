
process BAKTA_GETDB {

    tag "BAKTA DB v6.0 2025-02-24"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "$params.dbs", mode: 'copy'

    output:
    tuple path("bakta", type: "dir"), val("v6.0 2025-02-24"), emit: bakta_db

    script:
    """
    wget https://zenodo.org/record/14916843/files/db.tar.xz

    tar -xf db.tar.xz
    rm db.tar.xz
    mv db bakta

    echo "v6.0 2025-02-24" > bakta/VERSION.txt
    """
}
