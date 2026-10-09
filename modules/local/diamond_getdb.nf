process DIAMOND_GETDB {

    tag "DIAMOND ${db_name}"

    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/diamond:2.2.5--he361c42_0' :
        'biocontainers/diamond:2.2.5--he361c42_0' }"

    publishDir "${params.dbs}", mode: "copy"

    input:
    val db_name
    tuple path(protnlm_ids), val(protnlm_version)

    output:
    tuple path("${db_name}/db.filtered.dmnd"), env("VERSION"), emit: diamond_uniref_db
    path "${db_name}/VERSION.txt",                             emit: db_version_file
    path "${db_name}/README.md",                               emit: readme

    script:
    """
    wget "https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/${db_name}/RELEASE.metalink"

    mkdir -p ${db_name}

    grep -m1 "<version>" RELEASE.metalink | sed -e 's/.*<version>//' -e 's/<.*//' > ${db_name}/VERSION.txt

    VERSION=\$(cat ${db_name}/VERSION.txt)

    if [ "\$VERSION" != "${protnlm_version}" ]; then
        echo "UniRef release (\$VERSION) differs from TrEMBL release used for the ProtNLM list (${protnlm_version})" >&2
        exit 1
    fi

    wget "https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/${db_name}/${db_name}.fasta.gz"

    # Drop clusters whose representative name was predicted by Google ProtNLM (see README.md)
    # The cluster name is the representative's protein name, so checking RepID is enough
    gzip -dc ${db_name}.fasta.gz \\
        | awk 'NR==FNR{bad[\$1]; next}
               /^>/{n++; match(\$0, /RepID=[^ ]+/); keep = !(substr(\$0, RSTART+6, RLENGTH-6) in bad); if (!keep) d++}
               keep
               END{print n > "total.txt"; print d+0 > "removed.txt"}' ${protnlm_ids} - \\
        | diamond makedb --threads ${task.cpus} --db ${db_name}/db.filtered.dmnd

    cat <<-END_README > ${db_name}/README.md
    # ${db_name} DIAMOND database (ProtNLM-filtered)

    \\`db.filtered.dmnd\\` was built with \\`diamond makedb\\` from \\`${db_name}.fasta.gz\\`, UniProt release \$VERSION
    (https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/${db_name}/).

    ## Removed clusters

    Clusters whose representative's protein name was predicted by Google ProtNLM were removed.
    Those names are automatic predictions (evidence \\`ECO:0008006|Google:ProtNLM\\`) and are not reliable
    enough to transfer to hypothetical proteins.

    - Clusters in ${db_name}.fasta.gz: \$(cat total.txt)
    - Clusters removed: \$(cat removed.txt)
    - TrEMBL entries with ProtNLM names: \$(wc -l < ${protnlm_ids} | tr -d " ")

    END_README

    rm RELEASE.metalink ${db_name}.fasta.gz total.txt removed.txt
    """
}
