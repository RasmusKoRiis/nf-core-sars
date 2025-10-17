process MINIMAP2_ALIGN {
  tag { "${meta.id}" }
  publishDir "results/aln", mode: 'copy', overwrite: true

  container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
  cpus 4
  memory '4 GB'

  input:
  tuple val(meta), path(reads)    // List[*.fastq.gz]
  path  reference

  output:
  tuple val(meta),
        path("${meta.id}.bam"),
        path("${meta.id}.bam.bai"),
        emit: minimap2

  script:
  """
  minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads.join(' ')} \
    | samtools sort -@ ${task.cpus} -o ${meta.id}.bam
  # force BAI index (SARS-CoV-2 is tiny, so BAI is fine)
  samtools index -b ${meta.id}.bam
  """
}

