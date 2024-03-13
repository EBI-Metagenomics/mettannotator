process INTEPRO_ENTRY_LIST_GETDB {

    tag "InterPro Entry List 94.0"

    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("ipr-entry-lists/", type: "dir"), val("94.0"), emit: interpro_entry_list

    script:
    """
    wget https://ftp.ebi.ac.uk/pub/databases/interpro/releases/94.0/entry.list

    wget https://ftp.ebi.ac.uk/pub/databases/interpro/releases/94.0/ParentChildTreeFile.txt

    mkdir -p ipr-entry-lists

    mv entry.list ipr-entry-lists/

    mv ParentChildTreeFile.txt ipr-entry-lists/

    echo '94.0' > ipr-entry-lists/VERSION.txt

    """
}
