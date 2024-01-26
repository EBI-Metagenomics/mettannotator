process UNIFIRE {

    tag "${meta.prefix}"

    container "dockerhub.ebi.ac.uk/uniprot-public/unifire:2023.4"
    containerOptions "--bind unirule:/volume"

    input:
    tuple val(meta), path(faa, stageAs: "unirule/*")

    output:
    tuple val(meta), path("unirule/predictions_arba.out")                   , emit: arba
    tuple val(meta), path("unirule/predictions_unirule.out")                , emit: unirule
    tuple val(meta), path("unirule/predictions_unirule-pirsr.out")          , emit: pirsr
    tuple val(meta), path("unirule/${meta.prefix}.combined_predictions.out"), emit: unifire_all
    path("versions.yml")                                                    , emit: versions

    script:
    // This tool will only work using containers ATM.
    // it needs a specific folder to be mounted in order to work
    // we are mounting unirule in this case
    """
    echo "Pre-processed UniRule input file"
    prepare_unirule_input.py -i ${faa} -t ${meta.taxid} -o unirule

    # This is the provided docker running script
    /opt/scripts/bin/unifire-workflow.sh

    make_combined_unirule_output.py -i unirule -o unirule/${meta.prefix}.combined_predictions.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UniFIRE: 2023.4
    END_VERSIONS
    """
}
