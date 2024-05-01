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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Designed for dryad 
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
WorkflowDryad.initialise(params, log)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QUAST     } from '../modules/nf-core/quast'
include { MASHTREE  } from '../modules/nf-core/mashtree'
include { PARSNP    } from '../modules/nf-core/parsnp'
include { IQTREE    } from '../modules/nf-core/iqtree'
include { SNPDISTS  } from '../modules/nf-core/snpdists'

workflow DRYAD {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    INPUT_CHECK.out.reads
        .collect()
        .set { ch_input_reads }

    //
    // SUBWORKFLOW: Alignment Free
    //
    if (params.alignment_based == 'false') and (params.fasta == 'null') {
        MASHTREE (
            ch_input_reads
        )
    }
    ch_versions = ch_versions.mix(MASHTREE.out.versions.first())

    //
    // SUBWORKFLOW: Alignment Based
    //
    if (params.alignment_based == 'true') and (params.fasta != 'null') {
        ch_reference_fasta = params.fasta
        PARSNP (
            ch_input_reads
        )
        ch_versions = ch_versions.mix(PARSNP.out.versions.first()) 
        IQTREE (
            PARSNP.out.phylogeny
        )
        ch_versions = ch_versions.mix(IQTREE.out.versions.first())
        SNPDISTS (
            PARSNP.out.mblocks
        )
        ch_versions = ch_versions.mix(SNPDISTS.out.versions.first())
    }
}

// if (params.alignment_based == 'true') {
//     include {ALIGNMENT_BASED} from  './workflows/alignment_based'
// } else if (params.alignment_based == 'false') {
//     include {ALIGNMENT_FREE} from './workflows/alignment_free'
// }

// if (params.phoenix == 'false') {
    
// }