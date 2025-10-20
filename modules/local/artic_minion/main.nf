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
    path  bed optional true      
    path  ref optional true      

  output:
  tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus
  tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: artic_vcf

  script:
  def useBedRef = (bed && ref)            // both provided
  def bedRefFlag = useBedRef ? "--bed ${bed} --ref ${ref}" : ""
  def schemeFlag = useBedRef ? "" : "--scheme-directory ${scheme_dir} --scheme-name ${scheme_name} --scheme-version ${scheme_version}"
  def modelFlag  = params.artic_model ? "--model ${params.artic_model}" : ""
  """
  set -euo pipefail

  MODELDIR="\$PWD/clair3_models"
  mkdir -p "\$MODELDIR"
  artic_get_models --model-dir "\$MODELDIR" || true

  echo "ARTIC route: \$([ -n "${bedRefFlag}" ] && echo BED/REF || echo SCHEME)"
  [ -z "${bedRefFlag}" ] && {
    echo "scheme_dir: ${scheme_dir}"
    echo "scheme_name: ${scheme_name}"
    echo "scheme_version: ${scheme_version}"
    ls -l "${scheme_dir}/${scheme_name}/${scheme_version}" || true
  }

  artic minion \
    --normalise ${params.artic_normalise} \
    --threads ${task.cpus} \
    ${schemeFlag} \
    ${bedRefFlag} \
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
