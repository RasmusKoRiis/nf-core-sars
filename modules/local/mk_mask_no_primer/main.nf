process MK_MASK_NO_PRIMER {
    tag { "${meta.id}" }
    publishDir "results/qc", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    container 'bash:5.2'
    input:
    tuple val(meta), path(lowcov_bed)
    val  mask_primer_ends
    output:
    tuple val(meta), path("${meta.id}.mask.bed")

    when:
    !mask_primer_ends

    script:
    """
    cp ${lowcov_bed} ${meta.id}.mask.bed || touch ${meta.id}.mask.bed
    """
}