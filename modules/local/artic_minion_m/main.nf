process ARTIC_MINION_M {
  tag { "${meta.id}" }
  label 'process_high'
  publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'quay.io/artic/fieldbioinformatics:1.8.5'
  cpus { params.artic_threads }

  input:
    tuple val(meta), path(gp_fastq)

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: artic_vcf

  script:
  def modelFlag = params.artic_model ? "--model ${params.artic_model}" : ""

  """
  set -euo pipefail

  # Hard-coded Midnight v3
  REQ_NAME="ncov-2019_midnight"
  REQ_VER="v3.0.0"

  # Put scheme under the workdir so ARTIC can use --scheme-directory
  SCHEMERoot="\$PWD/primer_schemes"
  TARGET_DIR="\$SCHEMERoot/\$REQ_NAME/\$REQ_VER"
  BED_BN="\$REQ_NAME.scheme.bed"
  REF_BN="\$REQ_NAME.reference.fasta"
  BED_PATH="\$TARGET_DIR/\$BED_BN"
  REF_PATH="\$TARGET_DIR/\$REF_BN"

  echo "=== Ensuring Midnight v3 scheme is present ==="
  echo "Target: \$TARGET_DIR"
  mkdir -p "\$TARGET_DIR"

  need_dl=false
  [ -s "\$BED_PATH" ] || need_dl=true
  [ -s "\$REF_PATH" ] || need_dl=true

  fetch() {
    url="\$1"; out="\$2"
    if command -v curl >/dev/null 2>&1; then
      curl -fsSL --retry 5 --retry-delay 2 --connect-timeout 10 -o "\$out" "\$url"
    else
      wget -qO "\$out" "\$url"
    fi
  }

  if \$need_dl; then
    echo "Downloading Midnight v3 from quick-lab/primerschemes…"
    baseurl="https://raw.githubusercontent.com/quick-lab/primerschemes/master/\$REQ_NAME/\$REQ_VER"
    fetch "\$baseurl/\$BED_BN" "\$BED_PATH"
    fetch "\$baseurl/\$REF_BN" "\$REF_PATH"
  fi

  # Sanity checks & preview
  [ -s "\$BED_PATH" ] || { echo "ERROR: BED missing/empty: \$BED_PATH"; exit 2; }
  [ -s "\$REF_PATH" ] || { echo "ERROR: REF missing/empty: \$REF_PATH"; exit 2; }
  echo "== BED head =="
  head -n 3 "\$BED_PATH" || true
  echo "== REF header =="
  grep -m1 '^>' "\$REF_PATH" || true

  # Models
  MODELDIR="\$PWD/clair3_models"
  mkdir -p "\$MODELDIR"
  artic_get_models --model-dir "\$MODELDIR" || true

  # ARTIC run (scheme mode, hard-coded Midnight v3)
  artic minion \\
    --normalise ${params.artic_normalise} \\
    --threads ${task.cpus} \\
    --scheme-directory "\$SCHEMERoot" \\
    --scheme-name "\$REQ_NAME" \\
    --scheme-version "\$REQ_VER" \\
    --model-dir "\$MODELDIR" \\
    ${modelFlag} \\
    --read-file ${gp_fastq} \\
    ${meta.id}

  VCF=\$(ls ${meta.id}/clair3/*.pass.vcf.gz)
  tabix -p vcf "\$VCF" || true
  CONS=\$(ls ${meta.id}/*.consensus.fasta)

  cp "\$VCF" ${meta.id}.pass.vcf.gz
  [ -f "\$VCF.tbi" ] && cp "\$VCF.tbi" ${meta.id}.pass.vcf.gz.tbi || true
  cp "\$CONS" ${meta.id}.consensus.fasta
  """
}
