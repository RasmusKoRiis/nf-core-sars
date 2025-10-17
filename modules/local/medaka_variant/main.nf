process MEDAKA_VARIANT {
  tag { "${meta.id}" }
  publishDir "results/medaka", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/medaka:2.1.1--01dc988f451b713d'
  cpus 4
  memory '6 GB'

  input:
  tuple val(meta), path(bam)
  path  reference
  val   medaka_model   // e.g. r1041_e82_400bps_hac_variant_v4.2.0

  output:
  tuple val(meta), path("${meta.id}.vcf.gz"), emit: medaka_var

  script:
  """
  # 1) Convert coordinate-sorted BAM -> FASTQ (OK for haploid calling)
  samtools fastq ${bam} > ${meta.id}.tmp.fq

  # 2) Call variants (use -s to realign; better around indels)
  medaka_variant \
    -i ${meta.id}.tmp.fq \
    -r ${reference} \
    -o medaka_${meta.id} \
    -m ${medaka_model} \
    -s \
    -t ${task.cpus}

  # 3) Compress + index the correct VCF path for medaka >=2.x
  VCF="medaka_${meta.id}/medaka.annotated.vcf"
  bgzip -c "$VCF" > ${meta.id}.vcf.gz
  tabix -p vcf ${meta.id}.vcf.gz

  rm -f ${meta.id}.tmp.fq
  """
}
