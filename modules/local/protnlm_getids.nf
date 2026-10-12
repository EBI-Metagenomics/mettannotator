process PROTNLM_GETIDS {

    tag "UniProtKB/TrEMBL ${version} ProtNLM"

    label 'process_single'
    label 'process_long'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/diamond:2.2.5--he361c42_0' :
        'biocontainers/diamond:2.2.5--he361c42_0' }"

    input:
    val version

    output:
    tuple path("protnlm_ids.txt"), val(version), emit: protnlm_ids

    script:
    """
    BASE="https://ftp.ebi.ac.uk/pub/databases/uniprot"

    release_version() {
        wget -qO- "\$1/RELEASE.metalink" | grep -m1 "<version>" | sed -e 's/.*<version>//' -e 's/<.*//'
    }

    # Check the correct release ${version} and set the URL and unpacking command accordingly
    # Older releases keep a knowledgebase<version>.tar.gz file with the uniprot_trembl.dat.gz inside
    if [ "\$(release_version \$BASE/current_release/knowledgebase/complete)" = "${version}" ]; then
        URL="\$BASE/current_release/knowledgebase/complete/uniprot_trembl.dat.gz"
        UNPACK="cat"
    elif [ "\$(release_version \$BASE/previous_major_releases/release-${version}/knowledgebase)" = "${version}" ]; then
        URL="\$BASE/previous_major_releases/release-${version}/knowledgebase/knowledgebase${version}.tar.gz"
        UNPACK="tar -xzOf - uniprot_trembl.dat.gz"
    else
        echo "UniProt release ${version} not found in current_release nor previous_major_releases" >&2
        exit 1
    fi

    # Stream the TrEMBL file and extract IDs of entries predicted by Google ProtNLM
    wget -qO- "\$URL" \\
        | \$UNPACK \\
        | gzip -dc \\
        | awk '/^ID   /{id=\$2} /^DE   (RecName|SubName): Full=.*Google:ProtNLM/{print id}' > protnlm_ids.txt
    """
}
