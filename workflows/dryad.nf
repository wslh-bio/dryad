/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Figures out what params are nf-core and nextflow and parses them
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Checks to ensure input parameters exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) {if (param) { file(param, checkIfExists: true) } }

// Checks for mandatory parameters and puts it into a channel
if (params.input) {ch_input = file(params.input) } else { exit 1, 'Input samplesheet is not specified!'}

//Starting the workflow
WorkflowDryad.initialise(params, log)
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
include { ALIGNMENT_BASED   } from '../subworkflows/local/alignment_based'
include { ALIGNMENT_FREE    } from '../subworkflows/local/alignment_free'

workflow DRYAD {

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
    // Phoenix
    //
    if (!params.phoenix) {
        QUAST ( ch_input_reads )
        ch_versions = ch_versions.mix(QUAST.out.versions)
    }

    //
    // SUBWORKFLOW: Alignment Free
    //
    if (!params.alignment_based && !params.fasta) {
        ALIGNMENT_FREE (
            ch_input_reads
             )
    }

    //
    // SUBWORKFLOW: Alignment Based
    //
    if (params.alignment_based && params.fasta) {
        ALIGNMENT_BASED (
            ch_input_reads,
            params.fasta,
            params.outdir
            )
    }
}