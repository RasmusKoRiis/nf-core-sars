/*
Primer-only workflow for SARS-CoV-2 consensus FASTA files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This workflow ingests one or more multi-FASTA files, splits them so every
sequence becomes its own task, builds the primer database for the provided
scheme and runs the primer mismatch analysis only.
*/

include { PRIMER_MISMATCH } from '../modules/local/primer_mismatch/main'
include { BUILD_PRIMER_DB } from '../modules/local/build_primer_db/main'

params.fasta         = params.fasta ?: null
params.primer_bed    = params.primer_bed ?: null
params.primer_fasta  = params.primer_fasta ?: null
params.primer_set_name = params.primer_set_name ?: null

def parsePrimerFasta(fastaGlob) {
    Channel
        .fromPath(fastaGlob)
        .ifEmpty { error "No FASTA files found with pattern: ${fastaGlob}" }
        .map { file ->
            def sampleId = file.getBaseName().replaceAll(/\.(fa|fasta|fna)$/,'')
            def meta = [ id: sampleId ]
            tuple(meta, file)
        }
}

process SPLIT_PRIMER_FASTA {
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
    mkdir -p split && seqkit split2 --quiet -s 1 -O split "$fasta"
    for f in split/*.fasta; do
        base=\$(basename "\$f")
        mv "\$f" "${meta.id}_\$base"
    done
    """
}

workflow SARSCOVSEQPRIMERCHECK {

    main:

    def fastaGlob = params.fasta ?: error("Please provide --fasta pointing to the multi-FASTA input file(s)")
    def primerBedFile = params.primer_bed ? file(params.primer_bed) : error("Please provide --primer_bed for the primer scheme.")
    if( !primerBedFile.exists() ) {
        error "Primer BED file not found: ${primerBedFile}"
    }
    def primerFastaFile = params.primer_fasta ? file(params.primer_fasta) : error("Please provide --primer_fasta for the primer scheme.")
    if( !primerFastaFile.exists() ) {
        error "Primer FASTA file not found: ${primerFastaFile}"
    }
    def primerSetName = params.primer_set_name ?: primerBedFile.getBaseName().replaceAll(/\.bed$/,'')

    ch_fasta_raw = parsePrimerFasta(fastaGlob)

    SPLIT_PRIMER_FASTA(ch_fasta_raw)

    ch_single_seq = SPLIT_PRIMER_FASTA.out.split_fastas
        .flatMap { meta, files ->
            files.collect { f ->
                def prefix = meta.id ? "${meta.id}_" : ""
                def seq_id = f.getBaseName()
                    .replaceFirst("^${prefix}", '')
                    .replaceAll(/\.(fa|fasta|fna)$/,'')
                tuple([ id: seq_id ], f)
            }
        }

    BUILD_PRIMER_DB(
        Channel.value(primerBedFile),
        Channel.value(primerFastaFile),
        Channel.value(primerSetName)
    )

    PRIMER_MISMATCH(
        ch_single_seq,
        BUILD_PRIMER_DB.out.primer_db,
        Channel.value(params.runid)
    )
}
