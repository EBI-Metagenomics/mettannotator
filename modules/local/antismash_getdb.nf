process ANTISMASH_GETDB {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/antismash:7.1.0.1_2'

    output:
    path("antismash_db/"), emit: antismash_db

    script:
    """
    download-antismash-databases --database-dir antismash_db
    """
}
