process INTERPRO_ENTRY_LIST_GETDB {

    tag "InterPro Entry List ${params.interpro_entry_list_version}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("interpro_entry_list/", type: 'dir'), val(params.interpro_entry_list_version), emit: interpro_entry_list

    script:
    """
    wget https://ftp.ebi.ac.uk/pub/databases/interpro/releases/${params.interpro_entry_list_version}/entry.list

    wget https://ftp.ebi.ac.uk/pub/databases/interpro/releases/${params.interpro_entry_list_version}/ParentChildTreeFile.txt

    mkdir -p interpro_entry_list

    mv entry.list interpro_entry_list/

    mv ParentChildTreeFile.txt interpro_entry_list/

    echo '${params.interpro_entry_list_version}' > interpro_entry_list/VERSION.txt

    """
}
