
process DETECT_NCRNA {

    container 'quay.io/microbiome-informatics/genomes-pipeline.detect_rrna:v3'

    input:
    tuple val(meta), path(fasta)
    tuple path(rfam_ncrna_models), val(rfam_version)

    output:
    tuple val(meta), path('*.ncrna.deoverlap.tbl'), emit: ncrna_tblout
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmsearch: \$(cmsearch -h | grep -o '^# INFERNAL [0-9.]*' | sed 's/^# INFERNAL *//')
        Rfam version: $rfam_version
    END_VERSIONS
    """
}
