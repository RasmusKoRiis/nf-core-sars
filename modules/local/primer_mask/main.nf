process PRIMER_MASK {
    tag { "${meta.id}" }
    publishDir "results/qc", mode: 'copy', overwrite: true
    //errorStrategy 'ignore'
    
    container 'community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509'
    input:
    tuple val(meta), path(lowcov_bed)
    path primer_bed
    val  mask_primer_ends
    
    output:
    tuple val(meta), path("${meta.id}.mask.bed"), emit: primer_mask

    when:
    mask_primer_ends

    script:
    """
    # Merge lowcov + primer coords
    (cat ${lowcov_bed} ${primer_bed} 2>/dev/null || cat ${lowcov_bed}) \
      | sort -k1,1 -k2,2n | bedtools merge > ${meta.id}.mask.bed
    """
}