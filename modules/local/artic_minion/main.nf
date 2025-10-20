process ARTIC_MINION {
  tag { "${meta.id}" }
  label 'process_high'
  publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'quay.io/artic/fieldbioinformatics:1.6.0'
  cpus { params.artic_threads }

  input:
    tuple val(meta), path(gp_fastq)     // from guppyplex
    val   scheme_name
    val   scheme_version
    path  scheme_dir                    // can be a real directory or an empty placeholder

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus
    tuple val(meta), path("${meta.id}.pass.vcf.gz"),     path("${meta.id}.pass.vcf.gz.tbi"), emit: artic_vcf

  script:
  // If scheme_dir exists, we pass --scheme-directory. If not, ARTIC can auto-fetch.
  def schemeFlag = file(scheme_dir)?.exists() ? "--scheme-directory ${scheme_dir}" : ""
  def modelFlag  = params.artic_model ? "--model ${params.artic_model}" : ""

  """
  set -euo pipefail

  # Optional: fetch r10.4.1 models if needed (only once per env)
  # artic_get_models || true

  # Run the MinION pipeline (Clair3 is the only path now)
  artic minion \
    --normalise ${params.artic_normalise} \
    --threads ${task.cpus} \
    ${schemeFlag} \
    --scheme-name ${scheme_name} \
    --scheme-version ${scheme_version} \
    ${modelFlag} \
    --read-file ${gp_fastq} \
    ${meta.id}

  # Collect outputs (ARTIC writes into subfolders)
  VCF=\$(ls ${meta.id}/clair3/*.pass.vcf.gz)
  tabix -p vcf "\$VCF" || true
  CONS=\$(ls ${meta.id}/*.consensus.fasta)

  cp "\$VCF"                          ${meta.id}.pass.vcf.gz
  [ -f "\$VCF.tbi" ] && cp "\$VCF.tbi" ${meta.id}.pass.vcf.gz.tbi || true
  cp "\$CONS"                         ${meta.id}.consensus.fasta
  """
}
