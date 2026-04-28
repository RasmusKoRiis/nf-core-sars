# Offline Mode TODO

Working checklist for adding an explicit offline mode to `nf-core-sars`.

## Goal

Add an opt-in offline mode that keeps the current online behavior as the default. When offline mode is enabled, the wrapper and workflow must not require internet access during execution. All required pipeline code, plugins, containers, Nextclade datasets, ARTIC/Clair3 models, primer resources, lookup tables, and input files must already exist locally.

## Current External Dependencies

- Wrapper Git sync: `wrapper-sars-wgs-fixed.sh` runs `git pull` or `git clone` for `$HOME/ngs_scripts`.
- Wrapper pipeline sync: `wrapper-sars-wgs-fixed.sh` runs `nextflow pull RasmusKoRiis/nf-core-sars -r "$PIPELINE_BRANCH"`.
- Wrapper remote pipeline launch: `wrapper-sars-wgs-fixed.sh` runs `nextflow run RasmusKoRiis/nf-core-sars/main.nf -r "$PIPELINE_BRANCH"`.
- Nextflow plugin: `nextflow.config` declares `nf-schema@2.5.1`, which must be installed in the local Nextflow plugin cache before offline use.
- Nextflow custom config: `nextflow.config` references `params.custom_config_base = https://raw.githubusercontent.com/nf-core/configs/...`; this is already guarded by `NXF_OFFLINE`, but offline mode should set that explicitly.
- iGenomes config: `nextflow.config` defaults to `igenomes_ignore = false` and `igenomes_base = s3://ngi-igenomes/igenomes/`. This pipeline appears to use local SARS resources directly, so offline mode should set `--igenomes_ignore true`.
- Containers: modules use remote image names from Docker Hub, Quay, Seqera Wave/community library, and Galaxy Singularity URLs. Offline mode must rely on locally cached Docker images or local Apptainer/Singularity image files.
- Nextclade data: `modules/local/nextclade/main.nf` always runs `nextclade dataset get --name nextstrain/sars-cov-2/wuhan-hu-1/orfs`, which requires internet unless replaced with a local dataset path.
- ARTIC/Clair3 models: `modules/local/artic_minion_m/main.nf` runs `artic_get_models --model-dir "$MODELDIR" || true`; this tries to fetch models if they are not already available.
- Completion hooks: `subworkflows/local/pipeline_completion/main.nf` can call webhook/email paths when `--hook_url`, `--email`, or `--email_on_fail` are set. Offline wrapper mode should not pass these unless explicitly requested.
- Test profiles: `conf/test.config` and `conf/test_full.config` point to remote nf-core test data. Not relevant for production wrapper runs, but offline testing should use local test data.
- Wrapper SMB access: the wrapper uses `smbclient` to copy inputs and update dashboard DB files on the N-drive. This is local network storage, not internet, but it is still a network dependency. Decide whether "offline" means no internet only or no network at all.

## Implementation Tasks

- [ ] Define offline semantics.
  - Decision needed: does offline mode mean "no internet, local network/Samba still allowed", or "no network at all"?
  - Proposed default: offline means no internet. SMB/N-drive remains allowed because the wrapper currently depends on it for inputs and dashboard DBs.

- [x] Add pipeline parameters in `nextflow.config`.
  - Add `params.offline = false`.
  - Add `params.nextclade_dataset = null` for a pre-downloaded Nextclade dataset directory.
  - Add `params.artic_model_dir = null` for a pre-downloaded ARTIC/Clair3 model directory.
  - Consider `params.require_cached_containers = false` or document this as wrapper-only preflight instead.

- [x] Add the new parameters to `nextflow_schema.json`.
  - `offline`: boolean, default `false`.
  - `nextclade_dataset`: path to local Nextclade dataset directory, required when `offline=true` and workflows use `NEXTCLADE`.
  - `artic_model_dir`: path to local ARTIC/Clair3 model directory, required when `offline=true` and workflows use `ARTIC_MINION_M`.

- [x] Add pipeline-side config behavior for offline mode.
  - Done: when `params.offline` is true, remote nf-core custom config is not fetched.
  - Done: local custom config paths are still allowed in offline mode.
  - Done: iGenomes config is ignored in offline mode.
  - Remaining wrapper task: export `NXF_OFFLINE=true` before launch.
  - Decision: keep `nf-schema` validation available; the wrapper/preflight must verify the plugin is cached before offline use.

- [x] Update `modules/local/nextclade/main.nf`.
  - Done: online mode keeps current `nextclade dataset get`.
  - Done: offline mode requires `params.nextclade_dataset` and uses the staged dataset path as `--input-dataset`.
  - Done: missing or empty offline dataset paths fail with a clear error.
  - Done: outputs and `versions.yml` are unchanged.

- [x] Update `modules/local/artic_minion_m/main.nf`.
  - Done: online mode keeps current best-effort `artic_get_models`.
  - Done: offline mode requires `params.artic_model_dir` and uses the staged model directory as `MODELDIR` without running `artic_get_models`.
  - Done: missing or empty offline model paths fail with a clear error.
  - Deferred: strict `params.artic_model` layout validation, because ARTIC/Clair3 model directory layout should be confirmed on the server cache.
  - Done: current consensus, VCF, BAM, IUPAC, and `versions.yml` behavior is preserved.

- [x] Decide container offline strategy.
  - Done: first implementation uses Docker cached-image preflight in the wrapper, because the current wrapper uses `-profile docker,server`.
  - Done: `bin/prepare_offline_cache.sh` pulls the wrapper-required Docker images while online.
  - Done: wrapper offline mode requires all wrapper-required image tags to be available locally before starting Nextflow.
  - Future option: add Apptainer/Singularity cache support if the wrapper switches away from Docker.

- [ ] Pin unstable image tags.
  - Replace `latest` tags where possible, especially `docker.io/nextstrain/nextclade:latest` and `docker.io/cdcgov/irma:latest`.
  - This should be done as a separate tested change because it can alter tool versions and results.

- [x] Inventory required container images.
  - Done: wrapper-required images are listed in `wrapper-sars-wgs-fixed.sh` and `bin/prepare_offline_cache.sh`.
  - Done: optional images for other modules/workflow modes are available via `bin/prepare_offline_cache.sh --include-optional-images`.
  - Wrapper-required images:
    - `quay.io/nf-core/ubuntu:20.04`
    - `quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0`
    - `quay.io/artic/fieldbioinformatics:1.6.0`
    - `community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612`
    - `docker.io/rasmuskriis/nextclade-python`
    - `docker.io/nextstrain/nextclade:latest`
    - `docker.io/rasmuskriis/blast_python_pandas:amd64`
  - Apptainer/Singularity option: provide a local image cache directory and configure `singularity.cacheDir` or use process-level local image paths.
  - Optional non-wrapper images available through `--include-optional-images`:
    - `docker.io/cdcgov/irma:latest`
    - `community.wave.seqera.io/library/bcftools:1.22--a51ee80717c2467e`
    - `community.wave.seqera.io/library/bedtools_samtools:2932e857ecf6b5f2`
    - `community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c`
    - `community.wave.seqera.io/library/medaka:2.1.1--01dc988f451b713d`
    - `docker.io/library/bash:5.2`
    - `quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0`
    - `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`
    - `quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0`

- [x] Add wrapper offline option in `wrapper-sars-wgs-fixed.sh`.
  - Done: added `-o` to enable offline mode.
  - Done: prints offline mode in the startup summary.
  - Done: skips `git pull`, `git clone`, and `nextflow pull` in offline mode.
  - Done: runs a local pipeline checkout instead of `RasmusKoRiis/nf-core-sars/main.nf -r "$PIPELINE_BRANCH"` in offline mode.
  - Done: local pipeline path is configurable with `PIPELINE_DIR`, defaulting to `$HOME/.nextflow/assets/RasmusKoRiis/nf-core-sars`.
  - Done: exports `NXF_OFFLINE=true` before the offline Nextflow launch.
  - Done: passes `--offline true`, `--igenomes_ignore true`, `--nextclade_dataset <local path>`, and `--artic_model_dir <local path>`.

- [x] Add wrapper preflight checks for offline mode.
  - Done: verifies the local pipeline checkout exists.
  - Done: verifies the samplesheet, sample directory, primer directory, BED, reference, and lookup tables exist.
  - Done: verifies `nf-schema@2.5.1` is already cached.
  - Done: verifies the local Nextclade dataset path exists and is non-empty.
  - Done: verifies the local ARTIC model directory exists and is non-empty.
  - Done: verifies required Docker images are present locally with `docker image inspect`.
  - Done: fails early with actionable messages instead of letting Nextflow fail deep into execution.

- [x] Add offline wrapper input overrides.
  - Done: `--local-fastq-dir <dir>` uses a local FASTQ samples directory instead of copying FASTQs from the N-drive.
  - Done: `--local-samplesheet <csv>` uses a local samplesheet instead of `/mnt/tempdata/fastq/<RUN>.csv`.
  - Done: `--local-fasta <path-or-glob>` switches the wrapper to `fasta-workflow`.
  - Done: FASTA workflow skips primer/ARTIC-specific wrapper validation and ARTIC model preflight.
  - Done: local FASTQ/FASTA input paths are not removed by wrapper cleanup.

- [x] Add repo asset fallback for offline resources.
  - Done: mutation lookup tables resolve from `/mnt/tempdata/sars_db/assets` first and `$PIPELINE_DIR/assets` second.
  - Done: primer directories resolve from `/mnt/tempdata/sars_db/assets/<PRIMER>` first and `$PIPELINE_DIR/assets/<PRIMER>` second.
  - Done: resolved lookup table paths are passed to Nextflow.
  - Done: documented `PIPELINE_ASSETS_DIR` override.

- [x] Decide how offline resources are installed or refreshed.
  - Done: added `bin/prepare_offline_cache.sh` as the online cache-preparation procedure.
  - Done: script downloads/caches the Nextclade dataset once while online.
  - Done: script downloads/caches ARTIC/Clair3 models once while online.
  - Done: script pulls all wrapper-required Docker images once while online.
  - Done: script runs `nextflow plugin install nf-schema@2.5.1`.
  - Done: script runs `nextflow pull RasmusKoRiis/nf-core-sars -r <branch>`.

- [x] Add local resource defaults for the server.
  - Done: standard cache base defaults to `/mnt/tempdata/sars_db/assets/offline`.
  - Done: Nextclade dataset defaults to `/mnt/tempdata/sars_db/assets/offline/nextclade/sars-cov-2-wuhan-hu-1-orfs`.
  - Done: ARTIC models default to `/mnt/tempdata/sars_db/assets/offline/artic_models`.
  - Done: pipeline checkout defaults to `$HOME/.nextflow/assets/RasmusKoRiis/nf-core-sars`.
  - Done: paths are configurable through environment variables and cache-prep script flags.

- [ ] Make the wrapper release/version behavior offline-safe.
  - Online mode: keep `-b <branch>` and `nextflow pull`.
  - Offline mode: treat `-b <branch>` as informational unless the local checkout can be verified at that branch/tag.
  - Consider logging the local Git commit if the local pipeline directory is a Git repository.
  - Ensure `--release_version` reflects the offline pipeline version actually used.
  - Done: final FASTQ and FASTA reports now include a `Version Control Metadata` column with pipeline/Nextflow/offline mode, Docker image tags, and upstream process tool versions.

- [ ] Add tests.
  - Add module-level tests for `NEXTCLADE` online/offline command generation if practical.
  - Add module-level tests for `ARTIC_MINION_M` online/offline command generation if practical.
  - Add a wrapper dry-run/preflight test path, or a shellcheck-style validation if wrapper dry-run is not available.
  - Add a local stub workflow test with `--offline true` and local fake resource paths.

- [x] Add documentation.
  - Done: added `docs/offline-cache.md`.
  - Done: documented online cache preparation.
  - Done: documented wrapper offline usage and path overrides.
  - Done: documented wrapper preflight checks.
  - Remaining clarification: decide and document whether SMB/N-drive access counts as "offline" or not.

- [ ] Verification checklist before calling offline mode done.
  - Run online mode unchanged.
  - Run offline wrapper preflight with internet disabled or blocked.
  - Run a full offline workflow with already cached containers and local datasets.
  - Confirm logs contain no `nextflow pull`, `git pull`, `git clone`, `nextclade dataset get`, or `artic_get_models` in offline mode.
  - Confirm outputs match online mode for the same input run, except for timestamps and environment metadata.

## Suggested Implementation Order

1. Done: add `params.offline`, `params.nextclade_dataset`, and `params.artic_model_dir` to config/schema.
2. Done: patch `NEXTCLADE` to use a local dataset in offline mode.
3. Done: patch `ARTIC_MINION_M` to use a local model directory in offline mode.
4. Done: add wrapper `-o` flag, local pipeline launch, `NXF_OFFLINE=true`, and offline parameter passing.
5. Done: add wrapper offline preflight checks.
6. Done: document all required container images and add a cache-preparation procedure.
7. Next: run online regression, then offline end-to-end validation.

## Implementation Log

- 2026-04-27: Added pipeline offline parameters and schema entries.
- 2026-04-27: Updated config to skip remote custom configs and iGenomes when `--offline true`.
- 2026-04-27: Updated `NEXTCLADE` to use a staged local Nextclade dataset in offline mode.
- 2026-04-27: Updated `ARTIC_MINION_M` to use a staged local ARTIC/Clair3 model directory in offline mode.
- 2026-04-27: Validation performed with FASTA stub runs in normal and offline modes. Direct ARTIC stub command generation was checked, but local execution is blocked by the existing 12-CPU `process_high` label on this 8-CPU host.
- 2026-04-27: Added wrapper `-o` offline mode, local pipeline launch, `NXF_OFFLINE=true`, offline resource parameter passing, and preflight checks for cached plugin/resources/Docker images.
- 2026-04-27: Added `bin/prepare_offline_cache.sh` and `docs/offline-cache.md` to populate/document the offline cache while online.
- 2026-04-27: Added wrapper local input overrides for offline/local FASTQ and FASTA runs.
- 2026-04-27: Added wrapper fallback from `/mnt/tempdata/sars_db/assets` to `$PIPELINE_DIR/assets` for lookup tables and primer resources.
- 2026-04-27: Added final report `Version Control Metadata` column for both FASTQ and FASTA workflows.
