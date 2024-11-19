process ANNOTATE_GFF {

    tag "${meta.prefix}"

    label 'process_nano'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta),
        file(gff),
        file(eggnog_annotations_tsv),
        file(ncrna_tsv),
        file(trna_gff),
        file(crisprcas_hq_gff),
        file(amrfinder_tsv),
        file(antismash_gff),
        file(gecco_gff),
        file(dbcan_gff),
        file(df_gff),
        file(ips_annotations_tsv),     // empty in fast mode
        file(sanntis_annotations_gff), // empty in fast mode
        file(arba),                    // empty in fast mode
        file(unirule),                 // empty in fast mode
        file(pirsr)                    // empty in fast mode
    tuple path(interpro_entry_list), val(db_version)

    output:
    tuple val(meta), path("*_annotations.gff"),                                   emit: annotated_gff
    tuple val(meta), path("*_annotations_with_descriptions.gff"), optional: true, emit: annotated_desc_gff
    path "versions.yml",                                                          emit: versions

    script:
    def crisprcas_flag = "";
    def antismash_flag = "";
    def amrfinder_flag = "";
    def gecco_flag = "";
    def dbcan_flag = "";
    def df_flag = "";
    def ips_flag = "";
    def hypothetical_ips_flags = "";
    def sanntis_flag = "";
    if ( crisprcas_hq_gff ) {
        crisprcas_flag = "-c ${crisprcas_hq_gff} ";
    }
    if ( amrfinder_tsv ) {
        amrfinder_flag = "-a ${amrfinder_tsv}"
    }
    if ( antismash_gff ) {
        antismash_flag = "--antismash ${antismash_gff}"
    }
    if ( gecco_gff ) {
        gecco_flag = "--gecco ${gecco_gff}"
    }
    if ( dbcan_gff ) {
        dbcan_flag = "--dbcan ${dbcan_gff}"
    }
    if ( df_gff ) {
        df_flag = "--defense-finder ${df_gff}"
    }
    if ( ips_annotations_tsv ) {
        ips_flag = "-i ${ips_annotations_tsv}"
        hypothetical_ips_flags = [
            "--ipr-entries ${interpro_entry_list}/entry.list",
            "--ipr-hierarchy ${interpro_entry_list}/ParentChildTreeFile.txt",
            "--ipr-output ${ips_annotations_tsv}",
        ].join(" ")
    }
    if ( sanntis_annotations_gff ) {
        sanntis_flag = "-s ${sanntis_annotations_gff} ";
    }
    def hypothetical_tmp_gff = "${meta.prefix}_temp_with_unifire.gff"
    if ( params.fast ) {
        hypothetical_tmp_gff = "${meta.prefix}_temp.gff"
    }
    """
    annotate_gff.py \
    -g ${gff} \
    -e ${eggnog_annotations_tsv} \
    -r ${ncrna_tsv} \
    -t ${trna_gff} \
    -o ${meta.prefix}_temp.gff \
    ${crisprcas_flag} ${sanntis_flag} ${amrfinder_flag} \
    ${antismash_flag} ${gecco_flag} ${dbcan_flag} ${df_flag} ${ips_flag}

    if [ "${params.fast}" == "false" ]; then
        process_unifire_output.py \\
        -g ${meta.prefix}_temp.gff \\
        -a ${arba} \\
        -u ${unirule} \\
        -p ${pirsr} \\
        -o ${meta.prefix}_temp_with_unifire.gff
    fi

    add_hypothetical_protein_descriptions.py \\
    --eggnog-output ${eggnog_annotations_tsv} \\
    -i ${hypothetical_tmp_gff} \\
    -o ${meta.prefix}_annotations.gff ${hypothetical_ips_flags}

    if [ "${params.fast}" == "false" ]; then
        add_interpro_descriptions.py \\
        --ipr-entries ${interpro_entry_list}/entry.list \\
        -i ${meta.prefix}_annotations.gff \\
        -o ${meta.prefix}_annotations_with_descriptions.gff
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.prefix}_annotations.gff
    if [ "${params.fast}" == "false" ]; then
        touch ${meta.prefix}_annotations_with_descriptions.gff
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
