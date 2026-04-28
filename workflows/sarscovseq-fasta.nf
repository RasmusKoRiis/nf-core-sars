/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CSV_CONVERSION } from '../modules/local/csv_conversion/main'
include { NEXTCLADE      } from '../modules/local/nextclade/main'
include { REPORTFASTA    } from '../modules/local/report-fasta/main'
include { SPLIT_FASTA    } from '../modules/local/split_fasta/main'
include { TABLELOOKUP    } from '../modules/local/tablelookup/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parseFastaInputs(fastaPattern) {
    Channel
        .fromPath(fastaPattern, checkIfExists: true)
        .ifEmpty { error "No FASTA files found with pattern: ${fastaPattern}" }
        .map { fasta ->
            def sampleId = fasta.getBaseName().replaceAll(/\.(fa|fasta|fna)$/, '')
            tuple([id: sampleId], fasta)
        }
}

def resolveRequiredFile(pathValue, label) {
    if (!pathValue) {
        error "${label} was not provided."
    }
    def resolved = file(pathValue)
    if (!resolved.exists()) {
        error "${label} not found: ${pathValue}"
    }
    return resolved
}

def resolveRequiredDirectory(pathValue, label) {
    if (!pathValue) {
        error "${label} was not provided."
    }
    def resolved = file(pathValue)
    if (!resolved.exists()) {
        error "${label} not found: ${pathValue}"
    }
    if (!resolved.toFile().isDirectory()) {
        error "${label} is not a directory: ${pathValue}"
    }
    return resolved
}

def toPathList(paths) {
    paths instanceof List ? paths : [paths]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SARSCOVSEQFASTA {

    main:
    def fastaPattern = params.fasta ?: error('Please provide --fasta pointing to one or more FASTA files.')
    def spikeTable = resolveRequiredFile(params.spike, 'Spike resistance lookup table')
    def rdrpTable = resolveRequiredFile(params.rdrp, 'RdRp resistance lookup table')
    def clproTable = resolveRequiredFile(params.clpro, '3CLpro resistance lookup table')
    def offlineMode = params.offline.toString().toBoolean()
    def nextcladeDataset = offlineMode ? resolveRequiredDirectory(params.nextclade_dataset, 'Nextclade dataset directory') : []
    def reportVersionControlMetadata = [
        "pipeline=${workflow.manifest.name ?: 'nf-core-sars'}",
        "pipeline_version=${workflow.manifest.version ?: 'unknown'}",
        "nextflow=${workflow.nextflow.version}",
        "offline=${offlineMode}",
        "container_engine=${workflow.containerEngine ?: 'unknown'}",
        "docker_images=docker.io/rasmuskriis/nextclade-python|docker.io/nextstrain/nextclade:latest|docker.io/rasmuskriis/blast_python_pandas:amd64"
    ].join('; ')

    ch_fasta_raw = parseFastaInputs(fastaPattern)

    SPLIT_FASTA(ch_fasta_raw)

    ch_single_seq = SPLIT_FASTA.out.split_fastas
        .flatMap { meta, splitFastas ->
            toPathList(splitFastas).collect { splitFasta ->
                tuple([id: splitFasta.getBaseName()], splitFasta)
            }
        }

    NEXTCLADE(
        ch_single_seq,
        Channel.value(nextcladeDataset)
    )

    CSV_CONVERSION(NEXTCLADE.out.nextclade_csv)

    TABLELOOKUP(
        CSV_CONVERSION.out.nextclade_mutations,
        Channel.value(spikeTable),
        Channel.value(rdrpTable),
        Channel.value(clproTable)
    )

    ch_report_tool_versions = SPLIT_FASTA.out.versions
        .mix(NEXTCLADE.out.versions)
        .mix(CSV_CONVERSION.out.versions)
        .mix(TABLELOOKUP.out.versions)

    REPORTFASTA(
        CSV_CONVERSION.out.nextclade_stats_report.collect(),
        CSV_CONVERSION.out.nextclade_mutations_report.collect(),
        TABLELOOKUP.out.resistance_mutations_report.collect(),
        params.runid,
        ch_report_tool_versions.collect(),
        reportVersionControlMetadata
    )
}
