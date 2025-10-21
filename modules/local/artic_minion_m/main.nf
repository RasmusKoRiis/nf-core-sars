process ARTIC_MINION_M {
  tag { "${meta.id}" }
  label 'process_high'
  errorStrategy 'ignore'
  //publishDir "results/artic/${meta.id}", mode: 'copy', overwrite: true
  container 'community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612'
  cpus { params.artic_threads }

  input:
    tuple val(meta), path(gp_fastq)
    path  bed
    path  ref

  output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: artic_consensus
    path("${meta.id}.consensus.fasta"), emit: artic_consensus_report

  script:
  // Build the bash script as a single-quoted string (no Groovy interpolation),
  // then replace placeholders with actual values.
  def MODELOPT = (params.artic_model ? "--model ${params.artic_model}" : "").trim()

  def cmd = '''
set -euo pipefail

echo "=== ARTIC_MINION_M: BED/REF mode (consensus-only) ==="
echo "PWD: $(pwd)"
echo "BED: __BED__"
echo "REF: __REF__"
[ -s "__BED__" ] || { echo "ERROR: BED missing/empty: __BED__"; exit 2; }
[ -s "__REF__" ] || { echo "ERROR: REF missing/empty: __REF__"; exit 2; }

echo "BED head:"; head -n 3 "__BED__" || true
echo "REF header:"; grep -m1 '^>' "__REF__" || true

# Normalize BED to 7 columns (chrom, start, end, primername, pool, strand, sequence)
BED7=bed.normalized.bed
awk -F'\t' '
  BEGIN { OFS="\t" }
  /^#/ { print; next }                          # keep headers/comments
  NF>=7 { print; next }                         # already has sequence
  NF==6 {
    len = $3-$2; if (len<1) len=1;
    seq=""; for(i=0;i<len;i++) seq=seq "A";     # dummy sequence
    print $1,$2,$3,$4,$5,$6,seq; next
  }
  NF==5 {
    s = "+"; if (toupper($4) ~ /RIGHT/) s="-";  # infer strand from name
    len = $3-$2; if (len<1) len=1;
    seq=""; for(i=0;i<len;i++) seq=seq "A";
    print $1,$2,$3,$4,$5,s,seq; next
  }
  {
    print "ERROR: Unsupported BED line with " NF " columns: " $0 > "/dev/stderr";
    exit 11
  }
' "__BED__" > "$BED7"

echo "Normalized BED preview:"; head -n 3 "$BED7" || true

# Clair3 models (best effort; harmless if already present)
MODELDIR="$PWD/clair3_models"
mkdir -p "$MODELDIR"
artic_get_models --model-dir "$MODELDIR" || true

# Run ARTIC (direct bed/ref)
artic minion \
  --normalise __NORMALISE__ \
  --threads __THREADS__ \
  --bed "$BED7" \
  --ref "__REF__" \
  --model-dir "$MODELDIR" \
  __MODELOPT__ \
  --read-file __GPFASTQ__ \
  __METAID__

# Locate consensus robustly (ARTIC 1.6.x vs 1.8.x layout)
CONS=""
if [ -f "__METAID__.consensus.fasta" ]; then
  CONS="__METAID__.consensus.fasta"
elif ls __METAID__/*.consensus.fasta >/dev/null 2>&1; then
  for f in __METAID__/*.consensus.fasta; do
    CONS="$f"; break
  done
else
  echo "ERROR: Consensus fasta not found in known locations." >&2
  ls -lah || true
  exit 4
fi

# Emit expected filename only if needed (avoid copying onto itself)
if [ "$CONS" != "__METAID__.consensus.fasta" ]; then
  cp "$CONS" "__METAID__.consensus.fasta"
fi
'''

  // Substitute placeholders safely
  cmd = cmd
    .replace('__BED__', bed.toString())
    .replace('__REF__', ref.toString())
    .replace('__METAID__', meta.id.toString())
    .replace('__GPFASTQ__', gp_fastq.toString())
    .replace('__THREADS__', task.cpus.toString())
    .replace('__NORMALISE__', params.artic_normalise.toString())
    .replace('__MODELOPT__', MODELOPT)

  return cmd
}
