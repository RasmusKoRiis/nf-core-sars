process MAKE_DEPTH_MASK {
    tag { "${meta.id}" }
    publishDir "results/qc", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    container 'community.wave.seqera.io/library/pip_samtools_tool:1112661e1f5d0aae'
    input:
    tuple val(meta), path(bam)
    val   min_depth

    output:
    tuple val(meta), path("${meta.id}.lowcov.bed"), emit: depth_mask

    script:
    """
    # 1) Per-base depth
    samtools depth -a ${bam} > ${meta.id}.depth.txt

    # 2) Sites with depth < ${min_depth} -> BED
    awk 'BEGIN{OFS="\\t"} { if (\$3 < ${min_depth}) print \$1, \$2-1, \$2 }' ${meta.id}.depth.txt \
      | bedtools merge > ${meta.id}.lowcov.bed || touch ${meta.id}.lowcov.bed
    """
}