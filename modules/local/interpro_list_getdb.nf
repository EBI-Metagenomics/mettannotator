process INTEPRO_ENTRY_LIST_GETDB {

    tag "InterPro Entry List 105.0"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("interpro_entry_list/", type: "dir"), val("105.0"), emit: interpro_entry_list

    script:
    """
    wget https://ftp.ebi.ac.uk/pub/databases/interpro/releases/105.0/entry.list

    wget https://ftp.ebi.ac.uk/pub/databases/interpro/releases/105.0/ParentChildTreeFile.txt

    mkdir -p interpro_entry_list

    mv entry.list interpro_entry_list/

    mv ParentChildTreeFile.txt interpro_entry_list/

    echo '105.0' > interpro_entry_list/VERSION.txt

    """
}
