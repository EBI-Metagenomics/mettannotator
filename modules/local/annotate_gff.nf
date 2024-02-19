process ANNOTATE_GFF {

    tag "${meta.prefix}"

    label 'process_nano'

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta),
        file(gff),
        file(ips_annotations_tsv),
        file(eggnog_annotations_tsv),
        file(sanntis_annotations_gff),
        file(ncrna_tsv),
        file(crisprcas_hq_gff),
        file(amrfinder_tsv),
        file(antismash_gff),
        file(gecco_gff),
        file(dbcan_gff),
        file(df_gff),
        file(arba),
        file(unirule),
        file(pirsr)
    tuple path(interpro_entry_list), val(db_version)

    output:
    tuple val(meta), path("*_annotations.gff"), emit: annotated_gff
    path "versions.yml", emit: versions

    // For the version, I'm using the latest stable the genomes-annotation pipeline
    script:
    def eggnog_annotations_flag = ""
    def sanntis_flag = "";
    def crisprcas_flag = "";
    def amrfinder_flag = "";
    def antismash_flag = "";
    def gecco_flag = "";
    def dbcan_flag = "";
    def df_flag = "";
    if ( eggnog_annotations_tsv ) {
        eggnog_annotations_flag = "-e ${eggnog_annotations_tsv} "
    }
    if ( sanntis_annotations_gff ) {
        sanntis_flag = "-s ${sanntis_annotations_gff} ";
    }
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
    """
    annotate_gff.py \
    -g ${gff} \
    -i ${ips_annotations_tsv} \
    -r ${ncrna_tsv} \
    -o ${meta.prefix}_temp.gff \
    ${eggnog_annotations_flag} ${crisprcas_flag} ${sanntis_flag} ${amrfinder_flag} \
    ${antismash_flag} ${gecco_flag} ${dbcan_flag} ${df_flag}

    process_unifire_output.py \\
    -g ${meta.prefix}_temp.gff \\
    -a ${arba} \\
    -u ${unirule} \\
    -p ${pirsr} \\
    -o ${meta.prefix}_temp_with_unifire.gff

    add_hypothetical_protein_descriptions.py \\
    --ipr-entries ${interpro_entry_list}/entry.list \\
    --ipr-output ${ips_annotations_tsv} \\
    --eggnog-output ${eggnog_annotations_tsv} \\
    -i ${meta.prefix}_temp_with_unifire.gff \\
    -o ${meta.prefix}_annotations.gff


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.prefix}_annotated.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
