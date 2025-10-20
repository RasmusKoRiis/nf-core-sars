process ARTIC_MINION_M {
  tag { "${meta.id}" }
  label 'process_high'
  publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612'
  cpus { params.artic_threads }

  input:
    tuple val(meta), path(gp_fastq)
    path  bed
    path  ref

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus

  script:
  def modelFlag = params.artic_model ? "--model ${params.artic_model}" : ""
  """
  set -euo pipefail

  echo "=== ARTIC_MINION_M: BED/REF mode (consensus-only) ==="
  echo "PWD: \$(pwd)"
  echo "BED: ${bed}"
  echo "REF: ${ref}"
  [ -s "${bed}" ] || { echo "ERROR: BED missing/empty: ${bed}"; exit 2; }
  [ -s "${ref}" ] || { echo "ERROR: REF missing/empty: ${ref}"; exit 2; }

  echo "BED head:"; head -n 3 "${bed}" || true
  echo "REF header:"; grep -m1 '^>' "${ref}" || true

  # Normalize BED to 7 columns (chrom, start, end, primername, pool, strand, sequence)
  BED7=bed.normalized.bed
  awk -F'\\t' '
    BEGIN { OFS="\\t" }
    /^#/ { print; next }                          # keep headers/comments
    NF>=7 { print; next }                         # already has sequence
    NF==6 {
      len = \\$3-\\$2; if (len<1) len=1;
      seq=""; for(i=0;i<len;i++) seq=seq "A";     # dummy sequence
      print \\$1,\\$2,\\$3,\\$4,\\$5,\\$6,seq; next
    }
    NF==5 {
      s = "+"; if (toupper(\\$4) ~ /RIGHT/) s="-"; # infer strand from name
      len = \\$3-\\$2; if (len<1) len=1;
      seq=""; for(i=0;i<len;i++) seq=seq "A";
      print \\$1,\\$2,\\$3,\\$4,\\$5,s,seq; next
    }
    {
      print "ERROR: Unsupported BED line with " NF " columns: " \\$0 > "/dev/stderr";
      exit 11
    }
  ' "${bed}" > "\$BED7"

  echo "Normalized BED preview:"; head -n 3 "\$BED7" || true

  # Clair3 models (best effort; harmless if already present)
  MODELDIR="\$PWD/clair3_models"
  mkdir -p "\$MODELDIR"
  artic_get_models --model-dir "\$MODELDIR" || true

  # Run ARTIC (direct bed/ref). Even if we only want the consensus, ARTIC still
  # does its internal VCF steps; we just don’t emit them from the process.
  artic minion \\
    --normalise ${params.artic_normalise} \\
    --threads ${task.cpus} \\
    --bed "\$BED7" \\
    --ref "${ref}" \\
    --model-dir "\$MODELDIR" \\
    ${modelFlag} \\
    --read-file ${gp_fastq} \\
    ${meta.id}

  # Locate consensus robustly (ARTIC 1.6.x writes in CWD; 1.8.x under sample dir)
  CONS=""
  if [ -f "${meta.id}.consensus.fasta" ]; then
    CONS="${meta.id}.consensus.fasta"
  elif ls ${meta.id}/*.consensus.fasta >/dev/null 2>&1; then
    for f in ${meta.id}/*.consensus.fasta; do
      CONS="\\$f"; break
    done
  else
    echo "ERROR: Consensus fasta not found in known locations." >&2
    ls -lah || true
    exit 4
  fi

  # Emit expected filename only if needed (avoid copying onto itself)
  if [ "\$CONS" != "${meta.id}.consensus.fasta" ]; then
    cp "\$CONS" "${meta.id}.consensus.fasta"
  fi
  """
}
