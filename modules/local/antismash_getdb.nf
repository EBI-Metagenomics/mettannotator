process ANTISMASH_GETDB {

    tag "antiSMASH 8.0.1"

    publishDir "${params.dbs}", mode: 'copy'

    container "nf-core/antismash:8.0.1--pyhdfd78af_0"

    output:
    tuple path("antismash/", type: "dir"), val("8.0.1"), emit: antismash_db

    script:
    """
    download-antismash-databases --database-dir antismash

    echo "8.0.1" > antismash/VERSION.txt
    """
}
