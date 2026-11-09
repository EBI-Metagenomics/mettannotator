process DIAMOND_GETDB {

    tag "DIAMOND ${db_name}"

    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/diamond:2.2.5--he361c42_0' :
        'biocontainers/diamond:2.2.5--he361c42_0' }"

    publishDir "${params.dbs}", mode: "copy"

    input:
    val db_name

    output:
    tuple path("${db_name}/db.dmnd"), env("VERSION"), emit: diamond_uniref_db
    path "${db_name}/VERSION.txt",                    emit: db_version_file

    script:
    """
    wget "https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/${db_name}/${db_name}.fasta.gz"
    wget "https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/${db_name}/RELEASE.metalink"

    mkdir -p ${db_name}

    diamond makedb --threads ${task.cpus} --in ${db_name}.fasta.gz --db ${db_name}/db.dmnd

    grep -m1 "<version>" RELEASE.metalink | sed -e 's/.*<version>//' -e 's/<.*//' > ${db_name}/VERSION.txt

    VERSION=\$(cat ${db_name}/VERSION.txt)

    rm RELEASE.metalink
    """
}
