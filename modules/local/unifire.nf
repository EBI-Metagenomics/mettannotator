process UNIFIRE {

    tag "${meta.prefix}"

    container 'dockerhub.ebi.ac.uk/uniprot-public/unifire:2023.4'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("unirule/predictions_arba.out"), emit: arba
    tuple val(meta), path("unirule/predictions_unirule.out"), emit: unirule
    tuple val(meta), path("unirule/predictions_unirule-pirsr.out"), emit: pirsr

    script:
    """
    mkdir unirule && cp ${faa} unirule/${faa}
    echo "Copied input data for UniRule"
    python3 prepare_unirule_input.py -i unirule/${faa} -o unirule
    rm unirule/${faa}
    echo "Pre-processed UniRule input file"
    singularity run --bind unirule:/volume unifire:2023.4.sif
    echo "Finished running UniFire"
    """

}
