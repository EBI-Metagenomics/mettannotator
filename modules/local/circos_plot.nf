process CIRCOS_PLOT {
    tag "${meta.prefix}"

    label 'process_nano'

    container 'quay.io/microbiome-informatics/pycirclize:1.4.0'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.png"), emit: circos_plot, optional: true

    script:
    """
    circos_plot.py -p ${meta.prefix} -i ${gff} -o ${meta.prefix}_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyCirclize: \$(python -c "import pycirclize; print(pycirclize.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.prefix}_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyCirclize: \$(python -c "import pycirclize; print(pycirclize.__version__)")
    END_VERSIONS

    """
}
