/*
 * Predict tRNA genes
*/
process DETECT_TRNA {

    tag "${meta.prefix}"

    label 'process_nano'

    container 'quay.io/microbiome-informatics/genomes-pipeline.detect_rrna:v3.2'

    input:
    tuple val(meta), path(fasta), val(detected_kingdom)

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
    # and that causes issues as other detect_trna process will crash when the files are cleaned
    PROCESSTMP="\$PWD/tmp"
    mkdir -p "\$PWD/tmp"
    export TMPDIR="\${PROCESSTMP}"

    echo "[ Detecting tRNAs ]"
    tRNAscan-SE -${detected_kingdom} -Q \
    -m ${meta.prefix}_stats.out \
    -o ${meta.prefix}_trna.out \
    --gff ${meta.prefix}_trna_temp.gff \
    ${fasta}

    add_locus_tag_to_trna.py -i ${meta.prefix}_trna_temp.gff -o ${meta.prefix}_trna.gff
    
    parse_tRNA.py -i ${meta.prefix}_stats.out -o ${meta.prefix}_tRNA_20aa.out

    echo "Completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: \$(echo \$(tRNAscan-SE -h 2>&1) | grep -o "tRNAscan-SE [0-9].[0-9].[0-9]" | sed 's/tRNAscan-SE //')
    END_VERSIONS
    """
}
