/*
 * Predict bacterial 5S, 16S, 23S rRNA and tRNA genes
*/
process DETECT_RRNA {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.detect_rrna:v3.1'

    input:
    tuple val(meta), path(fasta)
    path cm_models

    output:
    path "results_folder/*.out", type: "file", emit: rrna_out_results
    path "results_folder/*.fasta", type: "file", emit: rrna_fasta_results
    path "results_folder/*.tblout.deoverlapped", emit: rrna_tblout_deoverlapped
    path "versions.yml", emit: versions


    // cmsearch_tblout_deoverlap version was taken from the container
    // it's using the same version that mgnify uses
    script:
    """
    shopt -s extglob

    RESULTS_FOLDER=results_folder
    FASTA=${fasta}
    CM_DB=${cm_models}

    FILENAME="${meta.prefix}"

    mkdir "\${RESULTS_FOLDER}"

    echo "[ Detecting rRNAs ] "

    for CM_FILE in "\${CM_DB}"/*.cm; do
        MODEL=\$(basename "\${CM_FILE}")
        echo "Running cmsearch for \${MODEL}..."
        cmsearch -Z 1000 \
            --hmmonly \
            --cut_ga --cpu ${task.cpus} \
            --noali \
            --tblout "\${RESULTS_FOLDER}/\${FILENAME}_\${MODEL}.tblout" \
            "\${CM_FILE}" "\${FASTA}" 1> "\${RESULTS_FOLDER}/\${FILENAME}_\${MODEL}.out"
    done

    echo "Concatenating results..."
    cat "\${RESULTS_FOLDER}/\${FILENAME}"_*.tblout > "\${RESULTS_FOLDER}/\${FILENAME}.tblout"

    echo "Removing overlaps..."
    cmsearch-deoverlap.pl \
    --maxkeep \
    --clanin "\${CM_DB}/ribo.claninfo" \
    "\${RESULTS_FOLDER}/\${FILENAME}.tblout"

    mv "\${FILENAME}.tblout.deoverlapped" "\${RESULTS_FOLDER}/\${FILENAME}.tblout.deoverlapped"

    echo "Parsing final results..."
    parse_rRNA-bacteria.py -i \
    "\${RESULTS_FOLDER}/\${FILENAME}.tblout.deoverlapped" 1> "\${RESULTS_FOLDER}/\${FILENAME}_rRNAs.out"

    rRNA2seq.py -d \
    "\${RESULTS_FOLDER}/\${FILENAME}.tblout.deoverlapped" \
    -i "\${FASTA}" 1> "\${RESULTS_FOLDER}/\${FILENAME}_rRNAs.fasta"

    echo "[ Detecting tRNAs ]"
    tRNAscan-SE -B -Q \
    -m "\${RESULTS_FOLDER}/\${FILENAME}_stats.out" \
    -o "\${RESULTS_FOLDER}/\${FILENAME}_trna.out" "\${FASTA}"

    parse_tRNA.py -i "\${RESULTS_FOLDER}/\${FILENAME}_stats.out" 1> "\${RESULTS_FOLDER}/\${FILENAME}_tRNA_20aa.out"

    echo "Completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: \$(echo \$(tRNAscan-SE -h 2>&1) | grep -o "tRNAscan-SE [0-9].[0-9].[0-9]" | sed 's/^# tRNAscan-SE //')
        cmsearch: \$(cmsearch -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//')
        cmsearch_tblout_deoverlap: 465b5026271
    END_VERSIONS
    """
}
