process DBCAN {

    tag "${meta.prefix}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/dbcan-5.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan-5.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(faa), path(gff)
    tuple path(dbcan_db, stageAs: "dbcan_db"), val(db_version)

    output:
    tuple val(meta), path("dbcan/substrate_prediction.tsv")  , emit: substrates
    tuple val(meta), path("dbcan/cgc_standard_out.tsv")      , emit: cgc
    tuple val(meta), path("dbcan/overview.tsv")              , emit: overview
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
        easy_substrate \\
        --threads ${task.cpus} \\
        --db_dir dbcan_db \\
        --output_dir dbcan \\
        --input_gff ${meta.prefix}_noseq.gff \\
        --input_raw_data ${faa} \\
        --gff_type prodigal \\
        --mode protein


    process_dbcan_result.py -i dbcan -o dbcan/${meta.prefix}_dbcan.gff -v 5.1.2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 5.1.2
        dbcan database: $db_version
    END_VERSIONS
    """
}
