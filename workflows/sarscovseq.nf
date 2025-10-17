
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

include { MINIMAP2_ALIGN              } from '../modules/local/minimap2_align/main'
include { MAKE_DEPTH_MASK             } from '../modules/local/make_depth_mask/main'
include { PRIMER_MASK                 } from '../modules/local/primer_mask/main'
include { MK_MASK_NO_PRIMER           } from '../modules/local/mk_mask_no_primer/main'
include { MEDAKA_VARIANT              } from '../modules/local/medaka_variant/main'
include { BCFTOOLS_CONSENSUS          } from '../modules/local/bcftools_consensus/main'







/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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
                log.warn "Skipping sample ${sampleId}: No FASTQ files found in ${files}"
                return null // Return null to allow filtering out later
            }

            // Creating a metadata map
            def meta = [ id: sampleId, single_end: true ]
            return tuple(meta, files)
        }
        .filter { it != null } // Remove null entries (samples with no FASTQ files)
}


workflow SARSCOVSEQ() {

    main:

    def currentDir = System.getProperty('user.dir')
    


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

    // MODULE MEDAK CONSENSUS 


    def min_depth = params.min_depth
    def mask_primer_ends   = params.mask_primer_ends
    def medaka_model   = params.medaka_model
 

    // Align
    MINIMAP2_ALIGN (
        CHOPPER.out.chopperfastq,
        Channel.value(file(params.reference))
    )

    // Build low-depth mask
    MAKE_DEPTH_MASK ( MINIMAP2_ALIGN.out.minimap2, min_depth )

    // Choose primer mask path or passthrough bed
    PRIMER_MASK (
    MAKE_DEPTH_MASK.out.lowcov,
    Channel.value(file(params.primer_bed)),
    mask_primer_ends
    )
    MK_MASK_NO_PRIMER (
    MAKE_DEPTH_MASK.out.lowcov,
    mask_primer_ends
    )
    def ch_mask = params.mask_primer_ends ? PRIMER_MASK.out.primer_mask : MK_MASK_NO_PRIMER.out.no_primer_mask

    // Medaka variants
    MEDAKA_VARIANT (
    MINIMAP2_ALIGN.out.minimap2,
    Channel.value(file(params.reference)),
    medaka_model
    )

    // --- Join VCF and mask by sample id (critical) ---
    def ch_vcf_tagged = MEDAKA_VARIANT.out.medaka_var.map { meta, vcf, tbi ->
    tuple(meta.id, meta, vcf, tbi)
    }
    def ch_mask_tagged = ch_mask.map { meta, bed ->
    tuple(meta.id, meta, bed)
    }

    // Join on id, then shape for bcftools
    def ch_consensus_in = ch_vcf_tagged.join(ch_mask_tagged).map { id, meta1, vcf, tbi, _id2, _meta2, bed ->
    tuple(meta1, vcf, tbi, bed)
    }

    // Consensus build
    BCFTOOLS_CONSENSUS (
    ch_consensus_in,
    Channel.value(file(params.reference))
    )


    //
    // MODULE: AMPLIGONE. SKIP IF MEDAKA KEEP IF IRMA
    //    

    AMPLIGONE (
        CHOPPER.out.chopperfastq, Channel.value(file(params.primerdir))
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
        //IRMA.out.bam
        MINIMAP2_ALIGN.out.minimap2
    )



    //
    // MODULE: NEXTCLADE
    //

    NEXTCLADE (
        //IRMA.out.amended_consensus
        BCFTOOLS_CONSENSUS.out.map { meta, fasta, report -> fasta } // pass consensus FASTA
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

    TABLELOOKUP (
        CSV_CONVERSION.out.nextclade_mutations, Channel.value(file(params.spike)), Channel.value(file(params.rdrp)), Channel.value(file(params.clpro))
    )

    //
    // MODULE: REPORT
    //

    def runid = params.runid
    def seq_instrument   = params.seq_instrument
    def primer   = "${params.primerdir}" 
    def release_version = params.release_version


    REPORT (
        CSV_CONVERSION.out.nextclade_stats_report.collect(), 
        CSV_CONVERSION.out.nextclade_mutations_report.collect(), 
        TABLELOOKUP.out.resistance_mutations_report.collect(),
        runid,
        release_version,
        seq_instrument,
        Channel.value(file(params.input)),
        primer,
        //IRMA.out.amended_consensus_report.collect()
        BCFTOOLS_CONSENSUS.out.map { meta, fasta, report -> report }.collect()
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
