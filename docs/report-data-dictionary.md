# Report data dictionary

This page documents the columns in the CSV report emitted by the `report` and
`report-fasta` modules. The report is a wide table: columns from the input
Nextclade and resistance files are carried through, so a run may contain fewer
gene columns than the complete list below. `NA` means unavailable or masked by
a quality gate; it is not a negative result.

| Column (or pattern) | Meaning and how it is made | Source used for the value | Important weakness / interpretation |
|---|---|---|---|
| `Sample` | Sample/sequence identifier; `seqName` or `SequenceID` is normalised to this join key. | Samplesheet, Nextclade `seqName`, resistance output | Duplicate identifiers are de-duplicated (first row retained), so inconsistent naming can merge records. |
| `RunID`, `Date`, `Instrument ID`, `Primer`, `Release Version` | Run metadata appended by the report process. Date is process date; other values are workflow parameters/assets. | Nextflow parameters and host clock | User-supplied values are not independently verified; date is not collection date. |
| `Barcode`, `PCR-PlatePosition`, `KonsCt` | Optional samplesheet identifiers and consensus Ct. | Input samplesheet | Manual/missing values are not independent QC; Ct is assay/platform dependent. |
| `clade`, `Nextclade_pango`, `partiallyAliased`, `clade_nextstrain`, `clade_who` | Nextclade clade/lineage classifications. | Nextclade output and reference dataset | Dataset-dependent; novel or recombinant sequences may be unresolved. |
| `coverage` | Fraction of reference genome covered. | Nextclade `coverage` | Genome-wide fraction hides local low-depth regions and base accuracy. |
| `cdsCoverage` | Per-gene coding-sequence coverage as `gene:fraction` pairs. | Nextclade `cdsCoverage` | Reference-coordinate and threshold dependent; missing genes can reflect parsing failures. |
| `NC_Genome_MixedSites`, `NC_Genome_QC`, `NC_Genome_frameShifts` | Nextclade mixed-site count, overall QC status, and frame-shift annotations. | Nextclade QC fields | Mixed sites have multiple possible causes; QC summary can conceal individual failures; indel calls are alignment dependent. |
| `<gene>_aaSubstitutions[_1..4]` | Amino-acid substitutions for SARS-CoV-2 genes; long S/ORF1a lists are split into numbered cells. | Nextclade fields grouped by `csv_conversion_nextclade.py` | Only reference-annotated calls are shown; low-coverage genes are masked. |
| `<gene>_aaDeletions`, `<gene>_aaInsertions` | Gene-specific amino-acid indels. | Nextclade fields grouped by `csv_conversion_nextclade.py` | Indel representation is alignment dependent; `No mutation` means none reported, not proof of absence. |
| `Spike_mAbs_inhibitors`, `3CLpro_inhibitors`, `RdRp_inhibitors` | Mutations matched against curated inhibitor/antibody tables. | `assets/*_inhibitors.csv` via `lookup_mutations.py` | Lookup coverage and literature may lag variants; unmatched is not proof of susceptibility. |
| `Spike_Fold`, `3CLpro_Fold`, `RdRp_Fold` | Maximum fold-change found for matched mutations. | Same inhibitor tables | Assay-specific and not automatically comparable; absent is not zero. |
| `DR_Res_Paxlovid`, `DR_Res_Remdesevir`, `DR_Res_mAbs` | Derived `Review`/`AANI` flag, or `NA` when unavailable. | Rule in `bin/report.py` | Triage only; does not model combinations, phenotype, dose, or treatment history. |
| `QC`, `NGS_QC_Sum` | Masking/low-coverage notes and flags such as `S:LC`, `Genome:MS`, `Genome:FS`. | `bin/report.py` and `bin/report_QC_calculation.py` | 80% gates and flag thresholds are policy choices; empty flags do not mean error-free. |
| `GISAID_Comment`, `GISAID_Kommentar` | `Review` when `NGS_QC_Sum` is non-empty; second is the same value with a German label. | `bin/report_QC_calculation.py` | Label duplication adds no evidence. |
| `NGS_Script_vers` | Semicolon-separated pipeline, tool, and Python version metadata. | Static metadata, `versions.yml`, Python runtime | May omit external database/laboratory changes. |

The report joins resistance summaries, Nextclade statistics/mutations, and
selected samplesheet fields by `Sample`, then applies 80% overall and 80% CDS
coverage masks before deriving resistance and QC flags. Retain the execution
report, trace, samplesheet, Nextclade dataset, primer asset, and inhibitor
tables with each report. Sources: [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/) and [GISAID](https://www.gisaid.org/).
