/*
 * Predict tRNA genes
*/
process DETECT_TRNA {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.detect_rrna:v3.2'

    input:
    tuple val(meta), path(fasta)
    tuple path(cm_models), val(rfam_version)

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

    # tRNAscan-SE needs a tmp folder otherwise it will use the base TMPDIR (with no subfolder)
    # and that causes issues as other detect_rrna process will crash when the files are cleaned
    PROCESSTMP="\$(mktemp -d)"
    export TMPDIR="\${PROCESSTMP}"
    # bash trap to clean the tmp directory
    trap 'rm -r -- "\${PROCESSTMP}"' EXIT


    echo "[ Detecting tRNAs ]"
    tRNAscan-SE -B -Q \
    -m "\${RESULTS_FOLDER}/\${FILENAME}_stats.out" \
    -o "\${RESULTS_FOLDER}/\${FILENAME}_trna.out" "\${FASTA}"

    parse_tRNA.py -i "\${RESULTS_FOLDER}/\${FILENAME}_stats.out" 1> "\${RESULTS_FOLDER}/\${FILENAME}_tRNA_20aa.out"

    echo "Completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: \$(echo \$(tRNAscan-SE -h 2>&1) | grep -o "tRNAscan-SE [0-9].[0-9].[0-9]" | sed 's/tRNAscan-SE //')
    END_VERSIONS
    """
}
