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
if (params.input) {ch_input = fo;e(params.input) } else { exit 1, 'Input samplesheet is not specified!'}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
if (params.phoenix == 'false') {
    include {QUAST} from './modules/nf-core/quast'
}
if (params.alignment_based == 'true') {
    include {ALIGNMENT_BASED} from  './workflows/alignment_based'
} else if (params.alignment_based == 'false') {
    include {ALIGNMENT_FREE} from './workflows/alignment_free'
}