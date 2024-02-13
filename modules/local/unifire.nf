process UNIFIRE {

    tag "${meta.prefix}"

    container "dockerhub.ebi.ac.uk/uniprot-public/unifire:2023.4"
    containerOptions "--bind unifire:/volume"

    input:
    tuple val(meta), path(faa, stageAs: "unifire/*")

    output:
    tuple val(meta), path("unifire/predictions_arba.out")                   , emit: arba
    tuple val(meta), path("unifire/predictions_unirule.out")                , emit: unirule
    tuple val(meta), path("unifire/predictions_unirule-pirsr.out")          , emit: pirsr
    tuple val(meta), path("unifire/${meta.prefix}.combined_predictions.out"), emit: unifire_all
    path("versions.yml")                                                    , emit: versions

    script:
    // This tool will only work using containers ATM.
    // it needs a specific folder to be mounted in order to work
    // we are mounting unifire in this case
    """
    echo "Pre-processed UniFIRE input file"
    prepare_unifire_input.py -i ${faa} -t ${meta.taxid} -o unifire

    # This is the provided docker running script
    /opt/scripts/bin/unifire-workflow.sh

    make_combined_unifire_output.py -i unifire -o unifire/${meta.prefix}.combined_predictions.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UniFIRE: 2023.4
    END_VERSIONS
    """
}
