/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Figures out what params are nf-core and nextflow and parses them
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Checks to ensure input parameters exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) {
    if (param) {
        file(param, checkIfExists: true)
        }
    }

// Checks for mandatory parameters and puts it into a channel
if (params.input) {
    ch_input = file(params.input) 
} 
else { 
    exit 1, 'Input samplesheet is not specified!'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Designed for dryad
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { QUAST             } from '../modules/local/quast'
include { QUAST_SUMMARY     } from '../modules/local/quast_summary'
include { ALIGNMENT_BASED   } from '../subworkflows/local/alignment_based'
include { ALIGNMENT_FREE    } from '../subworkflows/local/alignment_free'

workflow DRYAD {

    //
    // Error Handling
    //
    if (params.alignment_free && params.fasta && !params.alignment_based) {
        error("ERROR: An alignment free comparison does not use a reference fasta. Do you want to run an alignment based comparison instead?\nDryad terminating...")
        exit(1)
    }

    if (params.alignment_based && !params.fasta ) {
        error("ERROR: An alignment based comparison needs a reference fasta. If you want to run Dryad with a random reference picked by parsnp, use --fasta random .\nDryad terminating...")
        exit(1)
    }

    if (!params.alignment_based && !params.alignment_free) {
        error("ERROR: No alignment indicated. Please indicate which alignment to perform.")
        exit(1)
    }

    // Creating an empty channel to put version information into
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .set { ch_input_reads }

    // Adding version information
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // QC check for runs if skip quast
    //
    if (!params.skip_quast) {
        QUAST ( ch_input_reads )
        QUAST_SUMMARY (
            QUAST.out.transposed_report.collect()
            )
        ch_versions = ch_versions.mix(QUAST.out.versions)
    } else {
        quast_file = file("$baseDir/assets/empty.txt",checkIfExists:true)
    }

    //
    // Re-mapping channel to intake paths
    //
    ch_input_reads
        .map { sample, fasta ->
            fasta
        } // Produces queue channel of just fasta file paths in a list
        .collect()
        .set { ch_for_alignments }

    //
    // SUBWORKFLOW: Alignment Free
    //
    if (params.alignment_free) {
        ALIGNMENT_FREE (
            ch_for_alignments,
            params.task.cpus
             )
    }

    if (params.fasta == "random") {
        ch_fasta = file(params.random_file, checkIfExists:true)
    }
    else {
        ch_fasta = file(params.fasta, checkIfExists:true)
    }

    //
    // SUBWORKFLOW: Alignment Based
    //
    if (params.alignment_based && params.fasta) {
        if (!params.skip_quast) {
            ALIGNMENT_BASED (
                ch_for_alignments,
                ch_fasta,
                params.outdir,
                params.parsnp_partition,
                params.add_reference,
                INPUT_CHECK.out.csv,
                QUAST_SUMMARY.out.quast_tsv
                )
        }
        if (params.skip_quast) {
            ALIGNMENT_BASED (
                ch_for_alignments,
                ch_fasta,
                params.outdir,
                params.parsnp_partition,
                params.add_reference,
                INPUT_CHECK.out.csv,
                quast_file
                )
        }
    }
}