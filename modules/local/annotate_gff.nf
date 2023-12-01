process ANNOTATE_GFF {

    tag "${meta.prefix}"

    label 'process_light'

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta), file(gff), file(ips_annotations_tsv), file(eggnog_annotations_tsv), file(sanntis_annotations_gff),
    file(ncrna_tsv), file(crisprcas_hq_gff), file(amrfinder_tsv), path(arba), path(unirule), path(pirsr)

    output:
    tuple val(meta), path("*_annotated.gff"), emit: annotated_gff
    path "versions.yml", emit: versions

    // For the version, I'm using the latest stable the genomes-annotation pipeline
    script:
    def eggnog_annotations_flag = ""
    def sanntis_flag = "";
    def crisprcas_flag = "";
    def amrfinder_flag = "";
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
    """
    annotate_gff.py \
    -g ${gff} \
    -i ${ips_annotations_tsv} \
    -r ${ncrna_tsv} \
    -o ${meta.prefix}_temp.gff \
    ${eggnog_annotations_flag} ${crisprcas_flag} ${sanntis_flag} ${amrfinder_flag}

     process_unirule_output.py \
    -g ${meta.prefix}_temp.gff \
    -a ${arba} \
    -u ${unirule} \
    -p ${pirsr} \
    -o ${meta.prefix}_annotated.gff

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
