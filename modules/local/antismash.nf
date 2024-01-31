process ANTISMASH {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/antismash:7.1.0.1_2'

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${meta.prefix}_results/${meta.prefix}.gbk"), emit: gbk
    tuple val(meta), path("${meta.prefix}_antismash.tar.gz")          , emit: results_tarball
    tuple val(meta), path("${meta.prefix}_antismash.gff")             , emit: gff
    path "versions.yml"                                               , emit: versions

    script:
    """
    antismash \\
    -t bacteria \\
    -c ${task.cpus} \\
    --databases ${params.antismash_db} \\
    --output-basename ${meta.prefix} \\
    --genefinding-tool none \\
    --output-dir ${meta.prefix}_results \\
    ${gbk}

    tar -czf ${meta.prefix}_antismash.tar.gz ${meta.prefix}_results

    # To build the GFF3 file the scripts needs the regions.js file to be converted to json
    # In order to do that this process uses nodejs (using a patched version of the antismash container)

    echo ";var fs = require('fs'); fs.writeFileSync('./regions.json', JSON.stringify(recordData));" >> ${meta.prefix}_results/regions.js

    node ${meta.prefix}_results/regions.js

    antismash_to_gff.py \\
        -r regions.json -a \$(echo \$(antismash --version | sed 's/^antiSMASH //' )) \\
        -o ${meta.prefix}_antismash.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antiSMASH: \$(echo \$(antismash --version | sed 's/^antiSMASH //' ))
    END_VERSIONS
    """
}
