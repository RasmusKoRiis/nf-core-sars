# Offline Harddrive Installation

This guide creates a portable harddrive bundle for running the SARS workflow offline.

The offline computer must already have:

- Docker installed and running.
- Java and Nextflow available.
- The `NEXTFLOW` conda environment used by the wrapper, or an equivalent shell setup where `nextflow` works.

## 1. Prepare The Harddrive While Online

Mount the harddrive on a computer with internet access and run:

```bash
cd "$HOME/nf-core-sars"

bin/create_offline_harddrive_bundle.sh \
    --bundle-dir /media/$USER/SARS_OFFLINE \
    --branch master
```

Replace `/media/$USER/SARS_OFFLINE` with the actual mount path for the harddrive.

The script creates:

```text
SARS_OFFLINE/
├── assets/offline/
│   ├── artic_models/
│   └── nextclade/sars-cov-2-wuhan-hu-1-orfs/
├── docker/
│   ├── sars-docker-images.tar
│   └── sars-docker-images.txt
├── nextflow_home/
├── nf-core-sars/
├── README_OFFLINE_RUN.txt
└── run_offline_sars.sh
```

The script also downloads the Nextclade dataset, ARTIC/Clair3 models, `nf-schema@2.5.1`, the required Docker images, and exports the Docker images to `docker/sars-docker-images.tar`.

To also include optional images used by non-wrapper modules:

```bash
bin/create_offline_harddrive_bundle.sh \
    --bundle-dir /media/$USER/SARS_OFFLINE \
    --branch master \
    --include-optional-images
```

## 2. First Setup On The Offline Computer

Mount the harddrive and set a variable for the mount path:

```bash
export SARS_OFFLINE=/path/to/SARS_OFFLINE
```

Activate the environment if your setup requires it:

```bash
conda activate NEXTFLOW
```

If Conda is not installed under `~/miniconda3`, provide the Conda profile path when running:

```bash
export CONDA_PROFILE=/path/to/miniconda3/etc/profile.d/conda.sh
```

Load Docker images once:

```bash
docker load -i "$SARS_OFFLINE/docker/sars-docker-images.tar"
```

This imports the required images into the offline computer's local Docker cache.

## 3. Run FASTQ Workflow Offline

Use local FASTQ and samplesheet paths on the offline computer:

```bash
"$SARS_OFFLINE/run_offline_sars.sh" \
    -r <RUN> \
    -p <PRIMER> \
    -a sars \
    -y <YEAR> \
    --local-fastq-dir /path/to/fastq_pass \
    --local-samplesheet /path/to/samplesheet.csv
```

The FASTQ directory must contain barcode directories matching the `Barcode` column in the samplesheet.
By default, output is written to:

```text
$SARS_OFFLINE/results/<RUN>
```

By default, Nextflow work files are written to:

```text
$SARS_OFFLINE/work/<RUN>
```

To choose a different output directory:

```bash
"$SARS_OFFLINE/run_offline_sars.sh" \
    -r <RUN> \
    -p <PRIMER> \
    -a sars \
    -y <YEAR> \
    --local-fastq-dir /path/to/fastq_pass \
    --local-samplesheet /path/to/samplesheet.csv \
    --outdir /path/to/output/<RUN>
```

To choose a different Nextflow work directory:

```bash
"$SARS_OFFLINE/run_offline_sars.sh" \
    -r <RUN> \
    -p <PRIMER> \
    -a sars \
    -y <YEAR> \
    --local-fastq-dir /path/to/fastq_pass \
    --local-samplesheet /path/to/samplesheet.csv \
    --workdir /path/to/work/<RUN>
```

You can combine `--outdir` and `--workdir` when results and temporary work files should be stored on different disks.

Example layout:

```text
fastq_pass/
├── barcode01/
│   └── reads.fastq.gz
└── barcode02/
    └── reads.fastq.gz
```

Samplesheet columns must include at least:

```csv
SequenceID,Barcode
sample1,barcode01
sample2,barcode02
```

## 4. Run FASTA Workflow Offline

For already generated consensus FASTA files:

```bash
"$SARS_OFFLINE/run_offline_sars.sh" \
    -r <RUN> \
    -a sars \
    -y <YEAR> \
    --local-fasta "/path/to/*.fasta"
```

The FASTA workflow skips ARTIC/primer processing and runs downstream Nextclade/reporting only.

## 5. What The Launcher Sets

`run_offline_sars.sh` sets these paths automatically:

```bash
NXF_HOME="$SARS_OFFLINE/nextflow_home"
PIPELINE_DIR="$SARS_OFFLINE/nf-core-sars"
BASE_DIR="$SARS_OFFLINE/runtime"
TMP_DIR="$SARS_OFFLINE/runtime/fastq"
SARS_DATABASE="$SARS_OFFLINE/assets"
OFFLINE_BASE="$SARS_OFFLINE/assets/offline"
OFFLINE_OUTDIR_BASE="$SARS_OFFLINE/results"
OFFLINE_WORKDIR_BASE="$SARS_OFFLINE/work"
PIPELINE_ASSETS_DIR="$SARS_OFFLINE/nf-core-sars/assets"
```

Then it calls:

```bash
wrapper-sars-wgs-fixed.sh -o
```

## 6. Troubleshooting

If Docker images are missing:

```bash
docker load -i "$SARS_OFFLINE/docker/sars-docker-images.tar"
docker image ls
```

If `nf-schema@2.5.1` is missing, verify that the run uses the bundled Nextflow home:

```bash
export NXF_HOME="$SARS_OFFLINE/nextflow_home"
find "$NXF_HOME/plugins" -iname '*nf-schema*2.5.1*'
```

If Nextclade or ARTIC resources are missing:

```bash
find "$SARS_OFFLINE/assets/offline" -maxdepth 3 -type d
```

If the wrapper cannot find primer resources or mutation lookup tables, verify:

```bash
ls "$SARS_OFFLINE/nf-core-sars/assets"
```
