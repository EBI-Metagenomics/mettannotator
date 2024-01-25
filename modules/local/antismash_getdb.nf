process ANTISMASH_GETDB {

    tag "antiSMASH 7.1.0"

    publishDir "${params.dbs}", mode: 'copy'

    container 'quay.io/microbiome-informatics/antismash:7.1.0.1_2'

    output:
    tuple path("antismash/", type: "dir"), val("7.1.0"), emit: antismash_db

    script:
    """
    download-antismash-databases --database-dir antismash

    echo "7.1.0" > antismash/VERSION.txt
    """
}
