
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_sarscovseq_pipeline'

include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { AMPLIGONE                   } from '../modules/local/ampligone/main'
include { CHOPPER                     } from '../modules/local/chopper/main'
include { IRMA                        } from '../modules/local/irma/main'
include { NEXTCLADE                   } from '../modules/local/nextclade/main'
include { CSV_CONVERSION              } from '../modules/local/csv_conversion/main'
include { TABLELOOKUP                 } from '../modules/local/tablelookup/main'
include { REPORT                      } from '../modules/local/report/main'
include { DEPTH_ANALYSIS              } from '../modules/local/depth_analysis/main'







/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Function to parse the sample sheet
// Function to parse the sample sheet
def parseSampleSheet(sampleSheetPath) {
    return Channel
        .fromPath(sampleSheetPath)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            // Use 'SequenceID' and 'Barcode' based on the sample sheet
            if (!row.SequenceID || !row.Barcode) {
                error "Missing 'SequenceID' or 'Barcode' in the sample sheet row: ${row}"
            }

            // Use SequenceID as sample ID and Barcode to find files
            def sampleId = row.SequenceID
            def files = file("${params.samplesDir}/${row.Barcode}/*.fastq.gz")

            // Check if there are any files in the list
            if (!files || files.size() == 0) {
                error "No FastQ files for sample ${sampleId} found in ${files}"
            }

            // Creating a metadata map
            def meta = [ id: sampleId, single_end: true ]
            return tuple(meta, files)
        }
}


workflow SARSCOVSEQ() {

    main:

    def currentDir = System.getProperty('user.dir')
    def primerdir = "${currentDir}/${params.primerdir}"


    ch_sample_information = parseSampleSheet(params.input) // Use params.input directly
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
 

    ch_sample_information
        .map { meta, files ->
            tuple(meta, files.toList())
        }
    .set { read_input }


    //
    // MODULE: RUN CAT FASTQ
    //
    CAT_FASTQ (
        read_input
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: CHOPPER
    //  
    
    CHOPPER (
        CAT_FASTQ.out.reads
    )



    //
    // MODULE: AMPLIGONE.
    //    

    AMPLIGONE (
        CHOPPER.out.chopperfastq, primerdir
    )


    //
    // MODULE: IRMA
    //

    IRMA (
        AMPLIGONE.out.primertrimmedfastq
    )

    //
    // MODULE: DEPTH ANALYSIS
    //

    DEPTH_ANALYSIS (
        IRMA.out.bam
    )



    //
    // MODULE: NEXTCLADE
    //

    NEXTCLADE (
        IRMA.out.amended_consensus
    )

    //
    // MODULE: NEXTCLADE CONVERSION
    //

    CSV_CONVERSION (
        NEXTCLADE.out.nextclade_csv
    )

    //
    // MODULE: NEXTCLADE CONVERSION
    //

    def spike = "${currentDir}/${params.spike}"
    def rdrp = "${currentDir}/${params.rdrp}"
    def clpro = "${currentDir}/${params.clpro}"

    TABLELOOKUP (
        CSV_CONVERSION.out.nextclade_mutations, spike, rdrp, clpro
    )

    //
    // MODULE: REPORT
    //

    def runid = params.runid
    def seq_instrument   = params.seq_instrument  
    def samplesheet = "${currentDir}/assets/samplesheet.tsv"

    REPORT (
        CSV_CONVERSION.out.nextclade_stats_report.collect(), 
        CSV_CONVERSION.out.nextclade_mutations_report.collect(), 
        TABLELOOKUP.out.resistance_mutations_report.collect(),
        runid,
        seq_instrument,
        samplesheet
    )



    //
    // MODULE: Run FastQC
    //
    //FASTQC (
    //    read_input
    //)
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
