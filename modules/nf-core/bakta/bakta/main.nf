process BAKTA_BAKTA {
    tag "$meta.prefix"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.9.3--pyhdfd78af_0' :
        'biocontainers/bakta:1.9.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(detected_kingdom)
    tuple path(db), val(db_version)

    output:
    tuple val(meta), path("${prefix}.embl")             , emit: embl
    tuple val(meta), path("${prefix}.faa")              , emit: faa
    tuple val(meta), path("${prefix}.ffn")              , emit: ffn
    tuple val(meta), path("${prefix}.fna")              , emit: fna
    tuple val(meta), path("${prefix}.gbff")             , emit: gbk
    tuple val(meta), path("${prefix}.gff3")             , emit: gff
    tuple val(meta), path("${prefix}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}.txt")              , emit: txt
    tuple val(meta), path("${prefix}.svg")              , emit: svg
    tuple val(meta), path("${prefix}.png")              , emit: png
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.prefix}"
    """
    bakta \\
        $fasta \\
        $args \\
        --threads $task.cpus \\
        --prefix $prefix \\
        --locus-tag $prefix \\
        --db $db \\
        --keep-contig-headers \\
        --skip-trna \\
        --skip-tmrna \\
        --skip-rrna \\
        --skip-ncrna \\
        --skip-ncrna-region \\
        --skip-crispr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.embl
    touch ${prefix}.faa
    touch ${prefix}.ffn
    touch ${prefix}.fna
    touch ${prefix}.gbff
    touch ${prefix}.gff3
    touch ${prefix}.hypotheticals.tsv
    touch ${prefix}.hypotheticals.faa
    touch ${prefix}.tsv
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}
