/*
 * Predict tRNA genes
*/
process DETECT_TRNA {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.detect_rrna:v3.2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*_trna.out'), emit: trna_out
    tuple val(meta), path('*_stats.out'), emit: trna_stats
    tuple val(meta), path('*_tRNA_20aa.out'), emit: trna_count
    tuple val(meta), path('*_trna.gff'), emit: trna_gff
    path "versions.yml", emit: versions

    script:
    """
    shopt -s extglob

    # tRNAscan-SE needs a tmp folder otherwise it will use the base TMPDIR (with no subfolder)
    # and that causes issues as other detect_rrna process will crash when the files are cleaned
    PROCESSTMP="\$(mktemp -d)"
    export TMPDIR="\${PROCESSTMP}"
    # bash trap to clean the tmp directory
    trap 'rm -r -- "\${PROCESSTMP}"' EXIT


    echo "[ Detecting tRNAs ]"
    tRNAscan-SE -B -Q \
    -m ${meta.prefix}_stats.out \
    -o ${meta.prefix}_trna.out \
    --gff ${meta.prefix}_trna.gff \
    ${fasta}

    parse_tRNA.py -i ${meta.prefix}_stats.out -o ${meta.prefix}_tRNA_20aa.out

    echo "Completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: \$(echo \$(tRNAscan-SE -h 2>&1) | grep -o "tRNAscan-SE [0-9].[0-9].[0-9]" | sed 's/tRNAscan-SE //')
    END_VERSIONS
    """
}
