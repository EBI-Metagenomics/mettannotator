process DBCAN {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/dbcan:4.0.0'

    input:
    tuple val(meta), path(faa), path(gff),
    path dbcan_db

    output:
    tuple val(meta), path("dbcan/sub.prediction.out")        , emit: substrates
    tuple val(meta), path("dbcan/cgc_standard.out")          , emit: cgc
    tuple val(meta), path("dbcan/overview.txt")              , emit: overview
    tuple val(meta), path("dbcan/${meta.prefix}.dbcan.gff")  , emit: dbcan_gff
    path "versions.yml"                                      , emit: versions

    script:
    """
    while read line
    do
        if [[ \${line} == "##FASTA" ]]
        then
            break
        else
            echo "\$line"
        fi
    done < ${gff} > ${meta.prefix}_noseq.gff

    run_dbcan \\
        --dia_cpu ${task.cpus} \\
        --hmm_cpu ${task.cpus} \\
        --tf_cpu ${task.cpus} \\
        --db_dir ${dbcan_db} \\
        --out_dir dbcan \\
        --cgc_substrate \\
        --cluster ${meta.prefix}_noseq.gff \\
        ${faa} \\
        protein

    process_dbcan_result.py -i dbcan -o dbcan/${meta.prefix}.dbcan.gff -v 4.0.0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 4.0.0
    END_VERSIONS
    """
}
