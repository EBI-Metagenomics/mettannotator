process INTERPROSCAN_GETDB {

    tag "IPRS Scan ${params.interproscan_db_version}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("interproscan", type: "dir"), val("${params.interproscan_db_version}"), emit: interproscan_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${params.interproscan_db_version}/interproscan-${params.interproscan_db_version}-64-bit.tar.gz

    tar -xvzf interproscan-${params.interproscan_db_version}-64-bit.tar.gz

    mv interproscan-${params.interproscan_db_version}/ interproscan

    echo "${params.interproscan_db_version}" > interproscan/VERSION.txt

    rm interproscan-${params.interproscan_db_version}-64-bit.tar.gz
    """
}
