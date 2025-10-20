process ARTIC_MINION {
  tag { "${meta.id}" }
  label 'process_high'
  publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'quay.io/artic/fieldbioinformatics:1.8.5'
  cpus { params.artic_threads }

  input:
    tuple val(meta), path(gp_fastq)
    path  bed
    path  ref

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: artic_vcf

  script:
    def modelFlag = params.artic_model ? "--model ${params.artic_model}" : ""
  """
  set -euo pipefail

  echo "=== DEBUG: BED/REF mode only ==="
  echo "PWD: \$(pwd)"
  echo "BED path (staged): ${bed}"
  echo "REF path (staged): ${ref}"
  ls -l \$(dirname "${bed}") || true
  ls -l \$(dirname "${ref}") || true
  echo "BED head:"
  head -n 3 "${bed}" || true
  echo "REF header:"
  grep -m1 '^>' "${ref}" || true

  MODELDIR="\$PWD/clair3_models"
  mkdir -p "\$MODELDIR"
  artic_get_models --model-dir "\$MODELDIR" || true

  # sanity: make sure files exist and are non-empty
  [ -s "${bed}" ] || { echo "ERROR: BED missing/empty: ${bed}"; exit 2; }
  [ -s "${ref}" ] || { echo "ERROR: REF missing/empty: ${ref}"; exit 2; }

  artic minion \
    --normalise ${params.artic_normalise} \
    --threads ${task.cpus} \
    --bed "${bed}" \
    --ref "${ref}" \
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
