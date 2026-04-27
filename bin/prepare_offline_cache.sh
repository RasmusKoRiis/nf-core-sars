#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bin/prepare_offline_cache.sh [OPTIONS]

Populate the local cache required by wrapper-sars-wgs-fixed.sh -o.
Run this while online, on the same server/user that will run the pipeline offline.

Options:
  -h, --help                    Show this help message
  --base-dir <dir>              Offline cache base directory
                                Default: /mnt/tempdata/sars_db/assets/offline
  --nextclade-dir <dir>         Local Nextclade dataset directory
                                Default: <base-dir>/nextclade/sars-cov-2-wuhan-hu-1-orfs
  --artic-model-dir <dir>       Local ARTIC/Clair3 model directory
                                Default: <base-dir>/artic_models
  --pipeline-dir <dir>          Local pipeline checkout used by the wrapper
                                Default: $HOME/.nextflow/assets/RasmusKoRiis/nf-core-sars
  --pipeline-source <repo>      Nextflow pipeline source to pull
                                Default: RasmusKoRiis/nf-core-sars
  --branch <branch>             Pipeline branch/tag to cache
                                Default: master
  --nextflow <path>             Nextflow executable
                                Default: nextflow
  --plugin <id>                 Nextflow plugin to cache
                                Default: nf-schema@2.5.1
  --include-optional-images     Also pull images used by non-wrapper workflows/modules
  --skip-docker                 Do not pull Docker images
  --skip-nextclade              Do not download the Nextclade dataset
  --skip-artic-models           Do not download ARTIC/Clair3 models
  --skip-plugin                 Do not install/cache the Nextflow plugin
  --skip-pipeline               Do not run nextflow pull for the pipeline

Environment overrides:
  OFFLINE_BASE
  OFFLINE_NEXTCLADE_DATASET
  OFFLINE_ARTIC_MODEL_DIR
  PIPELINE_DIR
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

require_nonempty_dir() {
    local path="$1"
    local label="$2"
    [ -d "$path" ] || die "$label not found: $path"
    [ "$(find "$path" -mindepth 1 -maxdepth 1 | wc -l)" -gt 0 ] || die "$label is empty: $path"
}

OFFLINE_BASE="${OFFLINE_BASE:-/mnt/tempdata/sars_db/assets/offline}"
OFFLINE_NEXTCLADE_DATASET="${OFFLINE_NEXTCLADE_DATASET:-$OFFLINE_BASE/nextclade/sars-cov-2-wuhan-hu-1-orfs}"
OFFLINE_ARTIC_MODEL_DIR="${OFFLINE_ARTIC_MODEL_DIR:-$OFFLINE_BASE/artic_models}"
PIPELINE_DIR="${PIPELINE_DIR:-$HOME/.nextflow/assets/RasmusKoRiis/nf-core-sars}"
PIPELINE_SOURCE="RasmusKoRiis/nf-core-sars"
PIPELINE_BRANCH="master"
NEXTFLOW_BIN="nextflow"
PLUGIN_ID="nf-schema@2.5.1"
INCLUDE_OPTIONAL_IMAGES=false
DOCKER_ENABLED=true
NEXTCLADE_ENABLED=true
ARTIC_MODELS_ENABLED=true
PLUGIN_ENABLED=true
PIPELINE_ENABLED=true

NEXTCLADE_IMAGE="docker.io/nextstrain/nextclade:latest"
ARTIC_IMAGE="community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612"
NEXTCLADE_DATASET_NAME="nextstrain/sars-cov-2/wuhan-hu-1/orfs"

REQUIRED_DOCKER_IMAGES=(
    "quay.io/nf-core/ubuntu:20.04"
    "quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0"
    "quay.io/artic/fieldbioinformatics:1.6.0"
    "$ARTIC_IMAGE"
    "docker.io/rasmuskriis/nextclade-python"
    "$NEXTCLADE_IMAGE"
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
        --base-dir)
            OFFLINE_BASE="${2:?Missing value for --base-dir}"
            OFFLINE_NEXTCLADE_DATASET="$OFFLINE_BASE/nextclade/sars-cov-2-wuhan-hu-1-orfs"
            OFFLINE_ARTIC_MODEL_DIR="$OFFLINE_BASE/artic_models"
            shift 2
            ;;
        --nextclade-dir)
            OFFLINE_NEXTCLADE_DATASET="${2:?Missing value for --nextclade-dir}"
            shift 2
            ;;
        --artic-model-dir)
            OFFLINE_ARTIC_MODEL_DIR="${2:?Missing value for --artic-model-dir}"
            shift 2
            ;;
        --pipeline-dir)
            PIPELINE_DIR="${2:?Missing value for --pipeline-dir}"
            shift 2
            ;;
        --pipeline-source)
            PIPELINE_SOURCE="${2:?Missing value for --pipeline-source}"
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
        --plugin)
            PLUGIN_ID="${2:?Missing value for --plugin}"
            shift 2
            ;;
        --include-optional-images)
            INCLUDE_OPTIONAL_IMAGES=true
            shift
            ;;
        --skip-docker)
            DOCKER_ENABLED=false
            shift
            ;;
        --skip-nextclade)
            NEXTCLADE_ENABLED=false
            shift
            ;;
        --skip-artic-models)
            ARTIC_MODELS_ENABLED=false
            shift
            ;;
        --skip-plugin)
            PLUGIN_ENABLED=false
            shift
            ;;
        --skip-pipeline)
            PIPELINE_ENABLED=false
            shift
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

log "Offline cache settings"
echo "Offline base           : $OFFLINE_BASE"
echo "Nextclade dataset      : $OFFLINE_NEXTCLADE_DATASET"
echo "ARTIC model directory  : $OFFLINE_ARTIC_MODEL_DIR"
echo "Pipeline directory     : $PIPELINE_DIR"
echo "Pipeline source/branch : $PIPELINE_SOURCE / $PIPELINE_BRANCH"
echo "Nextflow executable    : $NEXTFLOW_BIN"
echo "Nextflow plugin        : $PLUGIN_ID"

if [ "$DOCKER_ENABLED" = true ] || [ "$NEXTCLADE_ENABLED" = true ] || [ "$ARTIC_MODELS_ENABLED" = true ]; then
    require_command docker
fi
if [ "$PLUGIN_ENABLED" = true ] || [ "$PIPELINE_ENABLED" = true ]; then
    require_command "$NEXTFLOW_BIN"
fi

if [ "$DOCKER_ENABLED" = true ]; then
    log "Pulling Docker images"
    docker_images=("${REQUIRED_DOCKER_IMAGES[@]}")
    if [ "$INCLUDE_OPTIONAL_IMAGES" = true ]; then
        docker_images+=("${OPTIONAL_DOCKER_IMAGES[@]}")
    fi
    for image in "${docker_images[@]}"; do
        echo "Pulling $image"
        docker pull "$image"
    done
fi

if [ "$NEXTCLADE_ENABLED" = true ]; then
    log "Downloading Nextclade dataset"
    mkdir -p "$(dirname "$OFFLINE_NEXTCLADE_DATASET")"
    docker run --rm \
        -v "$(dirname "$OFFLINE_NEXTCLADE_DATASET"):/offline-nextclade" \
        "$NEXTCLADE_IMAGE" \
        nextclade dataset get \
            --name "$NEXTCLADE_DATASET_NAME" \
            --output-dir "/offline-nextclade/$(basename "$OFFLINE_NEXTCLADE_DATASET")"
    require_nonempty_dir "$OFFLINE_NEXTCLADE_DATASET" "Nextclade dataset directory"
fi

if [ "$ARTIC_MODELS_ENABLED" = true ]; then
    log "Downloading ARTIC/Clair3 models"
    mkdir -p "$OFFLINE_ARTIC_MODEL_DIR"
    docker run --rm \
        -v "$OFFLINE_ARTIC_MODEL_DIR:/offline-artic-models" \
        "$ARTIC_IMAGE" \
        artic_get_models --model-dir /offline-artic-models
    require_nonempty_dir "$OFFLINE_ARTIC_MODEL_DIR" "ARTIC model directory"
fi

if [ "$PLUGIN_ENABLED" = true ]; then
    log "Caching Nextflow plugin"
    "$NEXTFLOW_BIN" plugin install "$PLUGIN_ID"
    if ! find "$HOME/.nextflow/plugins" -maxdepth 4 -iname "*${PLUGIN_ID%@*}*${PLUGIN_ID#*@}*" 2>/dev/null | grep -q .; then
        echo "WARNING: Could not verify $PLUGIN_ID in $HOME/.nextflow/plugins after install."
        echo "The wrapper offline preflight will still check for this before running."
    fi
fi

if [ "$PIPELINE_ENABLED" = true ]; then
    log "Caching Nextflow pipeline"
    "$NEXTFLOW_BIN" pull "$PIPELINE_SOURCE" -r "$PIPELINE_BRANCH"
    if [ -f "$PIPELINE_DIR/main.nf" ]; then
        echo "Verified local pipeline: $PIPELINE_DIR/main.nf"
    else
        echo "WARNING: Expected pipeline main.nf not found at $PIPELINE_DIR/main.nf"
        echo "If your local checkout is elsewhere, run the wrapper with PIPELINE_DIR=/path/to/nf-core-sars."
    fi
fi

log "Verifying required Docker image cache"
if [ "$DOCKER_ENABLED" = true ]; then
    missing_images=()
    for image in "${REQUIRED_DOCKER_IMAGES[@]}"; do
        if ! docker image inspect "$image" >/dev/null 2>&1; then
            missing_images+=("$image")
        fi
    done
    if [ "${#missing_images[@]}" -gt 0 ]; then
        printf 'Missing image: %s\n' "${missing_images[@]}" >&2
        die "Required Docker images are still missing."
    fi
fi

log "Offline cache preparation complete"
echo "Use with wrapper:"
echo "  PIPELINE_DIR=\"$PIPELINE_DIR\" \\"
echo "  OFFLINE_NEXTCLADE_DATASET=\"$OFFLINE_NEXTCLADE_DATASET\" \\"
echo "  OFFLINE_ARTIC_MODEL_DIR=\"$OFFLINE_ARTIC_MODEL_DIR\" \\"
echo "  ./wrapper-sars-wgs-fixed.sh -o -r <RUN> -p <PRIMER> -a sars -y <YEAR>"
