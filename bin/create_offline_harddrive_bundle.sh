#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bin/create_offline_harddrive_bundle.sh --bundle-dir <mounted-harddrive-dir> [OPTIONS]

Create a portable SARS offline bundle on a harddrive while this computer is online.
The offline computer is still expected to have Docker, Java/Nextflow, and the NEXTFLOW
conda environment installed.

Options:
  -h, --help                    Show this help message
  --bundle-dir <dir>            Harddrive directory to populate, for example /media/$USER/SARS_OFFLINE
  --branch <branch>             Pipeline branch/tag metadata to cache with Nextflow (default: master)
  --nextflow <path>             Nextflow executable to use while preparing cache (default: nextflow)
  --include-optional-images     Also export optional/non-wrapper module Docker images
  --skip-docker-save            Pull/cache images but do not create docker/sars-docker-images.tar
  --skip-cache-downloads        Only copy repo and write launch scripts; do not download datasets/models/plugins/images

Output layout:
  <bundle-dir>/nf-core-sars/       Pipeline checkout and scripts
  <bundle-dir>/assets/offline/     Nextclade dataset and ARTIC/Clair3 models
  <bundle-dir>/nextflow_home/      Portable Nextflow plugin/cache directory
  <bundle-dir>/docker/             Docker image tarball and image list
  <bundle-dir>/run_offline_sars.sh Convenience launcher for the offline computer
EOF
}

log() {
    printf '\n== %s ==\n' "$*"
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    local cmd="$1"
    command -v "$cmd" >/dev/null 2>&1 || die "Required command not found: $cmd"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

BUNDLE_DIR=""
PIPELINE_BRANCH="master"
NEXTFLOW_BIN="nextflow"
INCLUDE_OPTIONAL_IMAGES=false
DOCKER_SAVE_ENABLED=true
CACHE_DOWNLOADS_ENABLED=true

REQUIRED_DOCKER_IMAGES=(
    "quay.io/nf-core/ubuntu:20.04"
    "quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0"
    "quay.io/artic/fieldbioinformatics:1.6.0"
    "community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612"
    "docker.io/rasmuskriis/nextclade-python"
    "docker.io/nextstrain/nextclade:latest"
    "docker.io/rasmuskriis/blast_python_pandas:amd64"
)

OPTIONAL_DOCKER_IMAGES=(
    "docker.io/cdcgov/irma:latest"
    "community.wave.seqera.io/library/bcftools:1.22--a51ee80717c2467e"
    "community.wave.seqera.io/library/bedtools_samtools:2932e857ecf6b5f2"
    "community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c"
    "community.wave.seqera.io/library/medaka:2.1.1--01dc988f451b713d"
    "docker.io/library/bash:5.2"
    "quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0"
    "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
)

while [ "$#" -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --bundle-dir)
            BUNDLE_DIR="${2:?Missing value for --bundle-dir}"
            shift 2
            ;;
        --branch)
            PIPELINE_BRANCH="${2:?Missing value for --branch}"
            shift 2
            ;;
        --nextflow)
            NEXTFLOW_BIN="${2:?Missing value for --nextflow}"
            shift 2
            ;;
        --include-optional-images)
            INCLUDE_OPTIONAL_IMAGES=true
            shift
            ;;
        --skip-docker-save)
            DOCKER_SAVE_ENABLED=false
            shift
            ;;
        --skip-cache-downloads)
            CACHE_DOWNLOADS_ENABLED=false
            shift
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[ -n "$BUNDLE_DIR" ] || die "--bundle-dir is required."

if [ "$CACHE_DOWNLOADS_ENABLED" = true ]; then
    require_command docker
    require_command "$NEXTFLOW_BIN"
fi
if [ "$DOCKER_SAVE_ENABLED" = true ]; then
    require_command docker
fi
require_command rsync

BUNDLE_DIR="$(mkdir -p "$BUNDLE_DIR" && cd "$BUNDLE_DIR" && pwd)"
BUNDLE_PIPELINE_DIR="$BUNDLE_DIR/nf-core-sars"
BUNDLE_OFFLINE_BASE="$BUNDLE_DIR/assets/offline"
BUNDLE_NEXTFLOW_HOME="$BUNDLE_DIR/nextflow_home"
BUNDLE_DOCKER_DIR="$BUNDLE_DIR/docker"
BUNDLE_DOCKER_TAR="$BUNDLE_DOCKER_DIR/sars-docker-images.tar"

log "Bundle settings"
echo "Source repo       : $SOURCE_REPO_DIR"
echo "Bundle directory : $BUNDLE_DIR"
echo "Pipeline copy    : $BUNDLE_PIPELINE_DIR"
echo "Offline assets   : $BUNDLE_OFFLINE_BASE"
echo "Nextflow home    : $BUNDLE_NEXTFLOW_HOME"
echo "Docker tar       : $BUNDLE_DOCKER_TAR"
echo "Branch metadata  : $PIPELINE_BRANCH"

log "Copying pipeline repository"
mkdir -p "$BUNDLE_PIPELINE_DIR" "$BUNDLE_OFFLINE_BASE" "$BUNDLE_NEXTFLOW_HOME" "$BUNDLE_DOCKER_DIR"
rsync -a --delete \
    --exclude '.git/' \
    --exclude '.nextflow/' \
    --exclude 'work/' \
    --exclude 'results/' \
    --exclude 'out_sarsseq/' \
    "$SOURCE_REPO_DIR/" \
    "$BUNDLE_PIPELINE_DIR/"

chmod +x "$BUNDLE_PIPELINE_DIR/bin/prepare_offline_cache.sh" || true

if [ "$CACHE_DOWNLOADS_ENABLED" = true ]; then
    log "Downloading offline datasets, models, plugin, pipeline metadata, and Docker images"
    prepare_args=(
        --base-dir "$BUNDLE_OFFLINE_BASE"
        --pipeline-dir "$BUNDLE_PIPELINE_DIR"
        --branch "$PIPELINE_BRANCH"
        --nextflow "$NEXTFLOW_BIN"
    )
    if [ "$INCLUDE_OPTIONAL_IMAGES" = true ]; then
        prepare_args+=(--include-optional-images)
    fi
    NXF_HOME="$BUNDLE_NEXTFLOW_HOME" "$BUNDLE_PIPELINE_DIR/bin/prepare_offline_cache.sh" "${prepare_args[@]}"
else
    log "Skipping cache downloads"
fi

docker_images=("${REQUIRED_DOCKER_IMAGES[@]}")
if [ "$INCLUDE_OPTIONAL_IMAGES" = true ]; then
    docker_images+=("${OPTIONAL_DOCKER_IMAGES[@]}")
fi

printf '%s\n' "${docker_images[@]}" > "$BUNDLE_DOCKER_DIR/sars-docker-images.txt"

if [ "$DOCKER_SAVE_ENABLED" = true ]; then
    log "Exporting Docker images"
    docker save -o "$BUNDLE_DOCKER_TAR" "${docker_images[@]}"
fi

log "Writing offline launch helper"
cat > "$BUNDLE_DIR/run_offline_sars.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export NXF_HOME="${NXF_HOME:-$SCRIPT_DIR/nextflow_home}"
export PIPELINE_DIR="${PIPELINE_DIR:-$SCRIPT_DIR/nf-core-sars}"
export BASE_DIR="${BASE_DIR:-$SCRIPT_DIR/runtime}"
export TMP_DIR="${TMP_DIR:-$BASE_DIR/fastq}"
export SARS_DATABASE="${SARS_DATABASE:-$SCRIPT_DIR/assets}"
export OFFLINE_BASE="${OFFLINE_BASE:-$SCRIPT_DIR/assets/offline}"
export OFFLINE_OUTDIR_BASE="${OFFLINE_OUTDIR_BASE:-$SCRIPT_DIR/results}"
export PIPELINE_ASSETS_DIR="${PIPELINE_ASSETS_DIR:-$PIPELINE_DIR/assets}"

exec "$PIPELINE_DIR/wrapper-sars-wgs-fixed.sh" -o "$@"
EOF
chmod +x "$BUNDLE_DIR/run_offline_sars.sh"

log "Writing offline computer README"
cat > "$BUNDLE_DIR/README_OFFLINE_RUN.txt" <<EOF
SARS offline harddrive bundle

Prepared on: $(date -Iseconds)
Pipeline branch metadata: $PIPELINE_BRANCH

Before first offline run on a new computer:
  1. Mount this drive.
  2. Ensure Docker is running.
  3. Activate the NEXTFLOW conda environment if required:
       conda activate NEXTFLOW
     If Conda is not installed under ~/miniconda3, set:
       export CONDA_PROFILE=/path/to/miniconda3/etc/profile.d/conda.sh
  4. Load Docker images once:
       docker load -i "$BUNDLE_DOCKER_TAR"

Run FASTQ workflow:
  "$BUNDLE_DIR/run_offline_sars.sh" -r <RUN> -p <PRIMER> -a sars -y <YEAR> \\
      --local-fastq-dir /path/to/fastq_pass \\
      --local-samplesheet /path/to/samplesheet.csv

Output defaults to:
  "$BUNDLE_DIR/results/<RUN>"

To override:
  "$BUNDLE_DIR/run_offline_sars.sh" -r <RUN> -p <PRIMER> -a sars -y <YEAR> \\
      --local-fastq-dir /path/to/fastq_pass \\
      --local-samplesheet /path/to/samplesheet.csv \\
      --outdir /path/to/output/<RUN>

Run FASTA workflow:
  "$BUNDLE_DIR/run_offline_sars.sh" -r <RUN> -a sars -y <YEAR> \\
      --local-fasta "/path/to/*.fasta"

The launcher sets:
  NXF_HOME="$BUNDLE_NEXTFLOW_HOME"
  PIPELINE_DIR="$BUNDLE_PIPELINE_DIR"
  BASE_DIR="$BUNDLE_DIR/runtime"
  TMP_DIR="$BUNDLE_DIR/runtime/fastq"
  SARS_DATABASE="$BUNDLE_DIR/assets"
  OFFLINE_BASE="$BUNDLE_OFFLINE_BASE"
  OFFLINE_OUTDIR_BASE="$BUNDLE_DIR/results"
  PIPELINE_ASSETS_DIR="$BUNDLE_PIPELINE_DIR/assets"
EOF

log "Offline harddrive bundle complete"
echo "Next step on offline computer:"
echo "  docker load -i \"$BUNDLE_DOCKER_TAR\""
echo "  \"$BUNDLE_DIR/run_offline_sars.sh\" -r <RUN> -p <PRIMER> -a sars -y <YEAR> --local-fastq-dir <dir> --local-samplesheet <csv>"
