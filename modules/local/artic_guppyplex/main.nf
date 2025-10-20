process ARTIC_GUPPYPLEX {
  tag { "${meta.id}" }
  label 'process_medium'
  //publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'quay.io/artic/fieldbioinformatics:1.6.0'

  input:
    tuple val(meta), path(reads)    // reads is List[fastq.gz] for one sample

  output:
    tuple val(meta), path("${meta.id}_guppyplex.fastq"), emit: gp_fastq

  script:
  """
  set -euo pipefail
  # Concatenate all reads for this sample
  cat ${reads.join(' ')} > ${meta.id}.all.fq.gz

  # Length filtering (skip QC if only using 'pass' reads)
  artic guppyplex \
    --skip-quality-check \
    --min-length ${params.artic_min_len} \
    --max-length ${params.artic_max_len} \
    --prefix ${meta.id} \
    --directory .

  # guppyplex writes ${meta.id}_*.fastq; consolidate to a single file name
  cat ${meta.id}_*.fastq > ${meta.id}_guppyplex.fastq
  """
}
