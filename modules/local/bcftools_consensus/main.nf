process BCFTOOLS_CONSENSUS {
  tag { "${meta.id}" }
  publishDir "results/consensus", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/bcftools:1.22--a51ee80717c2467e'

  input:
  tuple val(meta), path(vcf), path(vcf_tbi), path(mask_bed)   // <-- accept .tbi and mask together
  path  reference

  output:
  tuple val(meta), path("${meta.id}.consensus.fasta"), path("${meta.id}.consensus.report.txt")

  script:
  """
  # Ensure VCF is indexed (belt-and-suspenders)
  if [ ! -f ${vcf}.tbi ]; then
    bcftools index -t ${vcf}
  fi

  # Apply variants onto reference
  bcftools consensus -f ${reference} ${vcf} > ${meta.id}.tmp.fa

  # Mask low-depth / primer regions with Ns
  if [ -s ${mask_bed} ]; then
    bcftools consensus -f ${meta.id}.tmp.fa -m ${mask_bed} ${vcf} > ${meta.id}.consensus.fasta
  else
    mv ${meta.id}.tmp.fa ${meta.id}.consensus.fasta
  fi

  printf "sample\\t${meta.id}\\nvcf\\t${vcf}\\nmask\\t${mask_bed}\\n" > ${meta.id}.consensus.report.txt
  """
}
