process MK_MASK_NO_PRIMER {
  tag { "${meta.id}" }
  publishDir "results/qc", mode: 'copy', overwrite: true

  input:
    tuple val(meta), path(lowcov_bed)
    val  mask_primer_ends

  when:
    !mask_primer_ends

  output:
    tuple val(meta), path("${meta.id}.mask.bed"), emit: no_primer_mask  // <-- name it

  container 'bash:5.2'
  script:
  """
  cp ${lowcov_bed} ${meta.id}.mask.bed || touch ${meta.id}.mask.bed
  """
}