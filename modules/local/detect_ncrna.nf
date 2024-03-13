
process DETECT_NCRNA {
    tag "${meta.prefix}"
    container 'quay.io/microbiome-informatics/genomes-pipeline.detect_rrna:v3.2'

    input:
    tuple val(meta), path(fasta)
    tuple path(rfam_ncrna_models), val(rfam_version)

    output:
    tuple val(meta), path('*.ncrna.deoverlap.tbl'), emit: ncrna_tblout
    tuple val(meta), path('*_rRNAs.out'), emit: rrna_out_results
    tuple val(meta), path('*_rRNAs.fasta'), emit: rrna_fasta_results
    path "versions.yml", emit: versions

    script:
    """
    cmscan \
    --cpu ${task.cpus} \
    --tblout overlapped_${fasta.baseName} \
    --hmmonly \
    --clanin ${rfam_ncrna_models}/Rfam.clanin \
    --fmt 2 \
    --cut_ga \
    --noali \
    -o /dev/null \
    ${rfam_ncrna_models}/Rfam.cm \
    ${fasta}

    # De-overlap #
    grep -v " = " overlapped_${fasta.baseName} > ${meta.prefix}.ncrna.deoverlap.tbl

    echo "Parsing final results..."
    parse_rRNA-bacteria.py \
    -s cmscan \
    -i ${meta.prefix}.ncrna.deoverlap.tbl \
    -o ${meta.prefix}_rRNAs.out

    rRNA2seq.py -d \
    ${meta.prefix}.ncrna.deoverlap.tbl \
    -s cmscan \
    -i ${fasta} \
    -o ${meta.prefix}_rRNAs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmscan: \$(cmscan -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//')
        Rfam version: $rfam_version
    END_VERSIONS
    """
}
