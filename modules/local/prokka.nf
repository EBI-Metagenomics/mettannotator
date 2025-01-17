process PROKKA {

    tag "${meta.prefix}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0' :
        'biocontainers/prokka:1.14.6--pl526_0' }"

    input:
    tuple val(meta), path(fasta), val(detected_kingdom)
    val(mode) // standard or compliant

    output:
    tuple val(meta), file("${meta.prefix}_prokka/${meta.prefix}.gff"), emit: gff
    tuple val(meta), file("${meta.prefix}_prokka/${meta.prefix}.faa"), emit: faa
    tuple val(meta), file("${meta.prefix}_prokka/${meta.prefix}.fna"), emit: fna
    tuple val(meta), file("${meta.prefix}_prokka/${meta.prefix}.gbk"), emit: gbk
    tuple val(meta), file("${meta.prefix}_prokka/${meta.prefix}.ffn"), emit: ffn
    tuple val(meta), file("${meta.prefix}_prokka/${meta.prefix}.txt"), emit: txt
    path "versions.yml" , emit: versions

    script:
    def compliant_flag = "";
    def rfam_flag = "";
    if ( mode == "compliant" ){
        compliant_flag = "--compliant"
        rfam_flag = "--rfam"
    }
    """
    # TMP folder issues in Prokka - https://github.com/tseemann/prokka/issues/402
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$PWD/tmp"

    # Disable the Java VM performane gathering tool, for improved performance
    export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

    cat ${fasta} | tr '-' ' ' > ${meta.prefix}_cleaned.fasta

    prokka ${meta.prefix}_cleaned.fasta \
    --cpus ${task.cpus} \
    --kingdom ${detected_kingdom} \
    --outdir ${meta.prefix}_prokka \
    --prefix ${meta.prefix} \
    --force \
    --locustag ${meta.prefix} \
    ${compliant_flag} \
    ${rfam_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
