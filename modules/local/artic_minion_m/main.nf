process ARTIC_MINION_M {
  tag { "${meta.id}" }
  label 'process_high'
  publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612'
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

  # === Pin Midnight v3 (Quick Lab primer schemes) ===
  REQ_NAME="ncov-2019_midnight"
  REQ_VER="v3.0.0"

  SCHEMERoot="\$PWD/primer_schemes"
  TARGET_DIR="\$SCHEMERoot/\$REQ_NAME/\$REQ_VER"
  BED_BN="\$REQ_NAME.scheme.bed"
  REF_BN="\$REQ_NAME.reference.fasta"
  BED_PATH="\$TARGET_DIR/\$BED_BN"
  REF_PATH="\$TARGET_DIR/\$REF_BN"

  echo "=== Ensuring scheme present: \$REQ_NAME/\$REQ_VER ==="
  echo "Target dir: \$TARGET_DIR"
  mkdir -p "\$TARGET_DIR"

  py_fetch() {
  python - "\$1" "\$2" <<'PY'
  import sys, time, ssl, urllib.request
  url, out = sys.argv[1], sys.argv[2]
  ctx = ssl.create_default_context()
  for i in range(5):               # simple retries
      try:
          with urllib.request.urlopen(url, context=ctx, timeout=30) as r:
              with open(out, "wb") as f:
                  f.write(r.read())
          sys.exit(0)
      except Exception as e:
          if i == 4:
              raise
          time.sleep(2 + i)
  PY
  }

  need_dl=false
  [ -s "\$BED_PATH" ] || need_dl=true
  [ -s "\$REF_PATH" ] || need_dl=true

  if \$need_dl; then
    echo "Downloading from quick-lab/primerschemes (raw GitHub)…"
    baseurl="https://raw.githubusercontent.com/quick-lab/primerschemes/master/\$REQ_NAME/\$REQ_VER"
    py_fetch "\$baseurl/\$BED_BN" "\$BED_PATH"
    py_fetch "\$baseurl/\$REF_BN" "\$REF_PATH"
  fi

  # Sanity checks
  [ -s "\$BED_PATH" ] || { echo "ERROR: BED missing/empty: \$BED_PATH"; exit 2; }
  [ -s "\$REF_PATH" ] || { echo "ERROR: REF missing/empty: \$REF_PATH"; exit 2; }
  echo "== BED head =="; head -n 3 "\$BED_PATH" || true
  echo "== REF header =="; grep -m1 '^>' "\$REF_PATH" || true

  # === Clair3 models (allowed to no-op if already present) ===
  MODELDIR="\$PWD/clair3_models"
  mkdir -p "\$MODELDIR"
  artic_get_models --model-dir "\$MODELDIR" || true

  # === ARTIC run (scheme mode, pinned Midnight v3) ===
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

  cp "\$VCF"                          ${meta.id}.pass.vcf.gz
  [ -f "\$VCF.tbi" ] && cp "\$VCF.tbi" ${meta.id}.pass.vcf.gz.tbi || true
  cp "\$CONS"                         ${meta.id}.consensus.fasta
  """


}
