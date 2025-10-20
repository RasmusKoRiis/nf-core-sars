process ARTIC_MINION {
  tag { "${meta.id}" }
  label 'process_high'
  publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'quay.io/artic/fieldbioinformatics:1.8.5'
  cpus { params.artic_threads }

  input:
  tuple val(meta), path(gp_fastq)
  val   scheme_name
  val   scheme_version
  path  scheme_dir

  output:
  tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus
  tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: artic_vcf

  // choose flags based on presence of params.bed/ref
  script:
  def haveBed   = params.bed && file(params.bed).exists()
  def haveRef   = params.ref && file(params.ref).exists()
  def bedFlag   = haveBed ? "--bed ${params.bed}" : ""
  def refFlag   = haveRef ? "--ref ${params.ref}" : ""
  def modelFlag = params.artic_model ? "--model ${params.artic_model}" : ""
  def schemeFlag = (!haveBed && !haveRef && file(scheme_dir)?.exists())
                   ? "--scheme-directory ${scheme_dir} --scheme-name ${scheme_name} --scheme-version ${scheme_version}"
                   : ""

  """
  set -euo pipefail

  MODELDIR="\$PWD/clair3_models"
  mkdir -p "\$MODELDIR"
  artic_get_models --model-dir "\$MODELDIR" || true

  artic minion \
    --normalise ${params.artic_normalise} \
    --threads ${task.cpus} \
    ${schemeFlag} \
    ${bedFlag} \
    ${refFlag} \
    --model-dir "\$MODELDIR" \
    ${modelFlag} \
    --read-file ${gp_fastq} \
    ${meta.id}

  VCF=\$(ls ${meta.id}/clair3/*.pass.vcf.gz)
  tabix -p vcf "\$VCF" || true
  CONS=\$(ls ${meta.id}/*.consensus.fasta)

  cp "\$VCF"                          ${meta.id}.pass.vcf.gz
  [ -f "\$VCF.tbi" ] && cp "\$VCF.tbi" ${meta.id}.pass.vcf.gz.tbi || true
  cp "\$CONS"                         ${meta.id}.consensus.fasta
  """
}
