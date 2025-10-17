// modules/local/bcftools_consensus/main.nf
process BCFTOOLS_CONSENSUS {
  tag { "${meta.id}" }
  publishDir "results/consensus", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/bcftools:1.22--a51ee80717c2467e'

  input:
    tuple val(meta), path(vcf), path(mask_bed)   
    path  reference

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), path("${meta.id}.consensus.report.txt")

  script:
  """
  # Ensure VCF is indexed (safe even if already indexed)
  bcftools index -t ${vcf}

  # Apply variants
  bcftools consensus -f ${reference} ${vcf} > ${meta.id}.tmp.fa

  # Mask Ns
  if [ -s ${mask_bed} ]; then
    bcftools consensus -f ${meta.id}.tmp.fa -m ${mask_bed} ${vcf} > ${meta.id}.consensus.fasta
  else
    mv ${meta.id}.tmp.fa ${meta.id}.consensus.fasta
  fi

  printf "sample\\t${meta.id}\\nvcf\\t${vcf}\\nmask\\t${mask_bed}\\n" > ${meta.id}.consensus.report.txt
  """
}
