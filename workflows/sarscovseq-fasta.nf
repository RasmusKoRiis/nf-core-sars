/*
FASTA-based SARS‑CoV‑2 analysis pipeline — per‑sequence mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This workflow ingests one or more **multi‑record FASTA** files and then
splits them so that every **single sequence** is pushed through the shared
modules individually:

  1. **SPLIT_FASTA** (new, inline) — converts each multi‑FASTA to one file
     per record.
  2. **NEXTCLADE** — lineage / QC on *one* sequence at a time.
  3. **CSV_CONVERSION** — tidy stats and mutation lists.
  4. **TABLELOOKUP** — antiviral resistance mapping.
  5. **REPORT** — combined PDF/HTML across all sequences.

The split step guarantees module isolation, prevents directory‑name clashes
and improves parallelism.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NEXTCLADE      } from '../modules/local/nextclade/main'
include { CSV_CONVERSION } from '../modules/local/csv_conversion/main'
include { TABLELOOKUP    } from '../modules/local/tablelookup/main'
include { REPORTFASTA    } from '../modules/local/report-fasta/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS (override with `--param value`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta           = params.fasta ?: null                   // --fasta "data/*.fasta"
params.runid           = params.runid ?: 'UNKNOWN_RUN'          // --runid 2025‑08‑06
params.seq_instrument  = params.seq_instrument ?: 'UNKNOWN'     // --seq_instrument "MiSeq"
params.primerdir       = params.primerdir ?: "${baseDir}/primer/"
params.release_version = params.release_version ?: 'dev'

params.spike           = params.spike  ?: "${baseDir}/resources/spike.tsv"
params.rdrp            = params.rdrp   ?: "${baseDir}/resources/rdrp.tsv"
params.clpro           = params.clpro  ?: "${baseDir}/resources/clpro.tsv"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Wrap each FASTA file into a (meta, path) tuple so that it matches the
// SPLIT_FASTA input interface.

def parseFasta(fastaGlob) {
    Channel
        .fromPath(fastaGlob)
        .ifEmpty { error "No FASTA files found with pattern: ${fastaGlob}" }
        .map { file ->
            def sampleId = file.getBaseName().replaceAll(/\.(fa|fasta|fna)$/,'')
            def meta = [ id: sampleId ]
            tuple(meta, file)
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: SPLIT_FASTA  🔪
    Split multi‑record FASTA into individual‑record FASTA files.
    Each output file is renamed:  <sampleId>_<seqHeader>.fasta
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process SPLIT_FASTA {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir ?: 'results'}/split_fastas", mode: 'copy', optional: true

    container 'quay.io/biocontainers/seqkit:2.8.1--h9ee0642_0'

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("${meta.id}_*.fasta"), emit: split_fastas

    script:
    """
    # Split into tmp dir to avoid clashes, then rename with sample prefix
    mkdir -p split && seqkit split2 --quiet -s 1 -O split "$fasta"
    for f in split/*.fasta; do
        base=\$(basename "\$f")
        mv "\$f" "${meta.id}_\$base"
    done
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SARSCOVSEQFASTA {

     main:

    // ───── Ensure user provided a FASTA input ─────
    def fastaGlob = params.fasta ?: error("Please provide --fasta <glob|file> … mou5 cin4 (謀前)!")

    // ───── Build input channel ─────
    ch_fasta_raw = parseFasta(fastaGlob)

    // ───── Split each multi‑FASTA into per‑sequence FASTAs ─────
    SPLIT_FASTA(ch_fasta_raw)

    // ───── Flatten so each sequence gets its own tuple & unique meta.id ─────
    ch_single_seq = SPLIT_FASTA.out.split_fastas
        .flatMap { meta, files ->
            files.collect { f ->
                def seq_id = f.getBaseName().replaceAll(/\.(fa|fasta|fna)$/,'')
                tuple([ id: seq_id ], f)
            }
        }

    // 1️⃣  Nextclade analysis (one seq per task)
    NEXTCLADE(
        ch_single_seq
    )

    // 2️⃣  CSV → TSV + QC
    CSV_CONVERSION(
        NEXTCLADE.out.nextclade_csv
    )

    // 3️⃣  Antiviral resistance lookup
    TABLELOOKUP(
        CSV_CONVERSION.out.nextclade_mutations,
        Channel.value(file(params.spike)),
        Channel.value(file(params.rdrp)),
        Channel.value(file(params.clpro))
    )

    def runid = params.runid
    def seq_instrument   = params.seq_instrument
    def primer   = "${params.primerdir}" 
    def release_version = params.release_version

    // 4️⃣  Final report (no IRMA consensus file in this path)
    REPORTFASTA(
        CSV_CONVERSION.out.nextclade_stats_report.collect(),
        CSV_CONVERSION.out.nextclade_mutations_report.collect(),
        TABLELOOKUP.out.resistance_mutations_report.collect(),
        runid
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
