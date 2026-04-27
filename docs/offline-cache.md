# Offline Cache Preparation

Offline mode is implemented as an opt-in wrapper mode. Run the cache preparation while the server still has internet access, then run the wrapper with `-o` when internet access is unavailable.

## Prepare Cache

From the pipeline repository:

```bash
bin/prepare_offline_cache.sh
```

By default this populates:

```text
/mnt/tempdata/sars_db/assets/offline/
├── artic_models/
└── nextclade/
    └── sars-cov-2-wuhan-hu-1-orfs/
```

It also pulls the Docker images required by `wrapper-sars-wgs-fixed.sh -o`, installs/caches `nf-schema@2.5.1`, and runs `nextflow pull RasmusKoRiis/nf-core-sars -r master`.

## Custom Paths

Use custom paths if the cache is stored elsewhere:

```bash
bin/prepare_offline_cache.sh \
    --base-dir /mnt/tempdata/sars_db/assets/offline \
    --pipeline-dir "$HOME/nf-core-sars" \
    --branch master
```

The wrapper must use the same paths:

```bash
PIPELINE_DIR="$HOME/nf-core-sars" \
OFFLINE_NEXTCLADE_DATASET="/mnt/tempdata/sars_db/assets/offline/nextclade/sars-cov-2-wuhan-hu-1-orfs" \
OFFLINE_ARTIC_MODEL_DIR="/mnt/tempdata/sars_db/assets/offline/artic_models" \
./wrapper-sars-wgs-fixed.sh -o -r <RUN> -p <PRIMER> -a sars -y <YEAR>
```

## Local Input Overrides

Offline mode can use local input paths instead of copying FASTQs from the N-drive.

For local FASTQ input, provide the barcode/FASTQ directory and a matching samplesheet:

```bash
PIPELINE_DIR="$HOME/nf-core-sars" \
./wrapper-sars-wgs-fixed.sh \
    -o \
    -r <RUN> \
    -p <PRIMER> \
    -a sars \
    -y <YEAR> \
    --local-fastq-dir /path/to/fastq_pass \
    --local-samplesheet /path/to/samplesheet.csv
```

The local FASTQ directory is passed as `--samplesDir`, so it should contain barcode directories matching the `Barcode` column in the samplesheet.

For local FASTA input, use `--local-fasta`. This switches the wrapper to the FASTA workflow and skips primer/ARTIC-specific setup:

```bash
PIPELINE_DIR="$HOME/nf-core-sars" \
./wrapper-sars-wgs-fixed.sh \
    -o \
    -r <RUN> \
    -a sars \
    -y <YEAR> \
    --local-fasta "/path/to/*.fasta"
```

## Optional Images

The default image set covers the current SARS WGS wrapper path. To also cache images used by other modules or workflow modes:

```bash
bin/prepare_offline_cache.sh --include-optional-images
```

## What Offline Mode Checks

The wrapper preflight checks:

- local pipeline checkout
- samplesheet and FASTQ sample directory for FASTQ workflow runs
- local FASTA input for FASTA workflow runs
- primer directory, BED, reference FASTA, and lookup tables for FASTQ workflow runs
- local Nextclade dataset
- local ARTIC/Clair3 model directory for FASTQ workflow runs
- cached `nf-schema@2.5.1`
- required Docker images via `docker image inspect`

If one of these is missing, the wrapper exits before starting Nextflow.
