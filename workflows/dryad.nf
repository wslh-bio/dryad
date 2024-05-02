/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check to ensure input parameters exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) {if (param) { file(param, checkIfExists: true) } }

// Check for mandatory parameters
if (params.input) {ch_input = file(params.input) } else { exit 1, 'Input samplesheet is not specified!'}

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
include { QUAST             } from '../modules/nf-core/quast'
include { ALIGNMENT_BASED   } '../subworkflows/local/alignment_based'
include { ALIGNMENT_FREE    } '../subworkflows/local/alignment_free'

workflow DRYAD {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    INPUT_CHECK.out.reads
        .collect()
        .set { ch_input_reads }

    //
    // Phoenix
    //
    if (!params.phoenix) {
        QUAST ( ch_input_reads )
    }

    //
    // SUBWORKFLOW: Alignment Free
    //
    if (!params.alignment_based && !params.fasta) {
        ALIGNMENT_FREE ( ch_input_reads )
    }

    //
    // SUBWORKFLOW: Alignment Based
    //
    if (params.alignment_based && params.fasta) {
        ALIGNMENT_BASED ( ch_input_reads )
    }
}

// if (params.alignment_based == 'true') {
//     include {ALIGNMENT_BASED} from  './workflows/alignment_based'
// } else if (params.alignment_based == 'false') {
//     include {ALIGNMENT_FREE} from './workflows/alignment_free'
// }

// if (params.phoenix == 'false') {
    
// }