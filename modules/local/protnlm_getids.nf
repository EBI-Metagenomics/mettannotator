process PROTNLM_GETIDS {

    tag "UniProtKB/TrEMBL ProtNLM"

    label 'process_single'
    label 'process_long'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/diamond:2.2.5--he361c42_0' :
        'biocontainers/diamond:2.2.5--he361c42_0' }"

    // Kept as a record of what was removed from the UniRef DIAMOND databases
    publishDir "${params.dbs}/protnlm", mode: "copy"

    output:
    tuple path("protnlm_ids.txt"), env("VERSION"), emit: protnlm_ids
    path "VERSION.txt",                            emit: version_file

    script:
    // Always the current release, the only one UniProt keeps as plain files
    """
    URL="https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete"

    wget -qO- "\$URL/RELEASE.metalink" | grep -m1 "<version>" | sed -e 's/.*<version>//' -e 's/<.*//' > VERSION.txt
    VERSION=\$(cat VERSION.txt)

    # IDs of entries whose name was predicted by Google ProtNLM.
    # Only the first top-level name counts, as that is the name UniRef uses for the cluster.
    wget -qO- "\$URL/uniprot_trembl.dat.gz" \\
        | gzip -dc \\
        | awk '/^ID   /{id=\$2; named=0}
               /^DE   (RecName|SubName): Full=/ && !named++ && /Google:ProtNLM/{print id}' > protnlm_ids.txt

    [ -s protnlm_ids.txt ] || { echo "No ProtNLM entries found in TrEMBL \$VERSION" >&2; exit 1; }
    """
}
