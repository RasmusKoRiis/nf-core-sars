# Offline Cache Preparation

Offline mode is implemented as an opt-in wrapper mode. Run the cache preparation while the server still has internet access, then run the wrapper with `-o` when internet access is unavailable.

## Prepare Cache

From the pipeline repository:

```bash
bin/prepare_offline_cache.sh
```

For a portable harddrive installation, prefer:

```bash
bin/create_offline_harddrive_bundle.sh --bundle-dir /media/$USER/SARS_OFFLINE
```

See [Offline harddrive installation](offline-harddrive-install.md).

By default this populates:

```text
/mnt/tempdata/sars_db/assets/offline/
├── artic_models/
└── nextclade/
    └── sars-cov-2-wuhan-hu-1-orfs/
```

It also pulls the Docker images required by `wrapper-sars-wgs-fixed.sh -o`, installs/caches `nf-schema@2.5.1`, and runs `nextflow pull RasmusKoRiis/nf-core-sars -r master`.

When `--offline` is used, the pipeline loads `conf/offline.config` and caps process CPU requests at 8 CPUs. Online runs and normal server runs keep their existing CPU settings.

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

To store Nextflow work files somewhere specific, pass `--workdir`:

```bash
PIPELINE_DIR="$HOME/nf-core-sars" \
./wrapper-sars-wgs-fixed.sh \
    -o \
    -r <RUN> \
    -p <PRIMER> \
    -a sars \
    -y <YEAR> \
    --workdir /path/to/work/<RUN>
```

For offline harddrive bundles, the launcher defaults this to `$SARS_OFFLINE/work/<RUN>`. Override the default base with `OFFLINE_WORKDIR_BASE` or set a per-run path with `--workdir`.

## Repo Asset Fallback

The wrapper uses `/mnt/tempdata/sars_db/assets` first for primer resources and mutation lookup tables. If a file or primer folder is missing there, it falls back to local pipeline assets:

```text
$PIPELINE_DIR/assets, then the assets/ directory next to wrapper-sars-wgs-fixed.sh
```

This means offline mode can still use bundled repo resources such as:

```text
assets/Spike_mAbs_inhibitors.csv
assets/RdRP_inhibitors.csv
assets/3CLpro_inhibitors.csv
assets/V5.4.2/
assets/VMIDT.2/
assets/primer_schemes/ncov-2019_midnight/v3.0.0/
```

If your local repo assets are somewhere else, override:

```bash
PIPELINE_ASSETS_DIR="/path/to/nf-core-sars/assets"
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
- primer directory, BED, reference FASTA, and lookup tables from either `/mnt/tempdata/sars_db/assets` or `$PIPELINE_DIR/assets`
- local Nextclade dataset
- local ARTIC/Clair3 model directory for FASTQ workflow runs
- cached `nf-schema@2.5.1`
- required Docker images via `docker image inspect`

If one of these is missing, the wrapper exits before starting Nextflow.
