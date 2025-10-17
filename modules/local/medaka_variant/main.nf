process MEDAKA_VARIANT {
  tag { "${meta.id}" }
  publishDir "results/medaka", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/medaka:2.1.1--01dc988f451b713d'
  cpus 4
  memory '6 GB'

  input:
  tuple val(meta), path(bam), path(bai)
  path  reference
  val   medaka_model   // e.g. r1041_e82_400bps_hac_variant_v4.2.0

 output:
  tuple val(meta),
        path("${meta.id}.vcf.gz"),
        path("${meta.id}.vcf.gz.tbi"),
        emit: medaka_var

  script:
  """
  samtools fastq ${bam} > ${meta.id}.tmp.fq

  medaka_variant \
    -i ${meta.id}.tmp.fq \
    -r ${reference} \
    -o medaka_${meta.id} \
    -m ${medaka_model} \
    -s \
    -t ${task.cpus}

  bgzip -c medaka_${meta.id}/medaka.annotated.vcf > ${meta.id}.vcf.gz
  tabix -p vcf ${meta.id}.vcf.gz

  rm -f ${meta.id}.tmp.fq
  """
}
