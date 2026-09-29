# Primer checking

The primer-checker integration runs PCR and NGS checks for SARS-CoV-2 and is
enabled by default in the FASTQ and FASTA workflows. Disable it with
`--primer_check false`.

Checks run on the final `${runid}.fasta` emitted by the report process. Both
FASTQ and FASTA workflows publish this combined consensus file. PCR checks use
only the PCR schemes in a unified JSON; unrelated NGS panel assets in that JSON
are not required. NGS checks use the selected sequencing scheme separately.

The FASTQ and FASTA wrappers default to `/mnt/tempdata/sars_db/pcr-primers` for
the PCR database and pass it explicitly to Nextflow. Override it with
`-P /path/to/database` or the `PRIMER_CHECK_PCR` environment variable (`-P`
takes precedence).

For direct pipeline runs, supply `--primer_check_pcr /path/to/database`. NGS checks
use the selected primer scheme directory; override it with
`--primer_check_ngs_dir /path/to/scheme`. FASTA runs also need a scheme supplied
through `--primerdir` or `--primer_check_ngs_dir` for NGS checks.

See the [shared integration guide](https://github.com/RasmusKoRiis/primer-checker/blob/main/docs/PIPELINE_INTEGRATION.md)
for database layouts, wrapper overrides, latest-image deployment, ignored errors,
output interpretation and synthetic tests. The default container must be
published before first production use. Primer tasks use `errorStrategy 'ignore'`
and publish CSV/HTML under `primer_check/`; inspect `task_status.csv` for failures.

The local standalone harness is `tests/primer_check/main.nf`, with its own
`nextflow.config`. It can exercise the module using synthetic manifests without
running the sequencing pipeline or installing nf-core plugins. The canonical
harness and automated tests are maintained in primer-checker.
