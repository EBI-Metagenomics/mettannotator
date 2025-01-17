process DBCAN {

    tag "${meta.prefix}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/dbcan:4.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan:4.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(faa), path(gff)
    tuple path(dbcan_db, stageAs: "dbcan_db"), val(db_version)

    output:
    tuple val(meta), path("dbcan/substrate.out")             , emit: substrates
    tuple val(meta), path("dbcan/cgc_standard.out")          , emit: cgc
    tuple val(meta), path("dbcan/overview.txt")              , emit: overview
    tuple val(meta), path("dbcan/${meta.prefix}_dbcan.gff")  , emit: dbcan_gff
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
        --db_dir dbcan_db \\
        --out_dir dbcan \\
        --cgc_substrate \\
        --cluster ${meta.prefix}_noseq.gff \\
        ${faa} \\
        protein

    process_dbcan_result.py -i dbcan -o dbcan/${meta.prefix}_dbcan.gff -v 4.1.2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 4.1.2
        dbcan database: $db_version
    END_VERSIONS
    """
}
