// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// include { QUAST                  } from '../modules/nf-core/quast/main'
// include { ALIGNMENT_BASED        } from '../subworkflows/local/alignment_based'
// include { ALIGNMENT_BASED        } from '../subworkflows/local/alignment_'
// include { paramsSummaryMap       } from 'plugin/nf-validation'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dryad_pipeline'

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     RUN MAIN WORKFLOW
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// workflow DRYAD {

//     if (params.alignment_based == false) { 
//         ALIGNMENT_FREE (ch_samplesheet)
// }
//     else if (params.alignment_based == true) {
//         ALIGNMENT_BASED (ch_samplesheet)
//     }
// }
//     // take:
//     // ch_samplesheet // channel: samplesheet read in from --input

//     // main:

//     // ch_versions = Channel.empty()
//     // ch_multiqc_files = Channel.empty()

//     // //
//     // // MODULE: Run QUAST if samples are not from Phoenix
//     // //
//     // if (!params.phoenix) {
//     //     QUAST (
//     //         ch_samplesheet.sample
//     //     )
//     // }

//     // // Updating versions channel
//     // ch_versions = ch_versions.mix(QUAST.out.versions.first())

//     // //
//     // // MODULE: Run subworkflow of alignment based or alignment free?
//     // //

//     // // If alignment free workflow is desired
//     // if (!params.alignment_based) {
//     //     ALIGNMENT_FREE (
//     //         INPUT_CHECK.out.reads
//     //     )
//     // }

//     // // If alignment based workflow is desired
//     // else {
//     //     ALIGNMENT_BASED (
//     //         INPUT_CHECK.out.reads
//     //     )
//     // }toList

//     // //
//     // // Collate and save software versions
//     // //
//     // softwareVersionsToYAML(ch_versions)
//     //     .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
//     //     .set { ch_collated_versions }

//     // //
//     // // MODULE: MultiQC
//     // //
//     // ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//     // ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
//     // ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
//     // summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
//     // ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
//     // ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
//     // ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
//     // ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//     // ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
//     // ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

//     // MULTIQC (
//     //     ch_multiqc_files.collect( ),
//     //     ch_multiqc_config.toList(),
//     //     ch_multiqc_custom_config.toList(),
//     //     ch_multiqc_logo.toList()
//     // )

//     // emit:
//     //     multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
//     //     versions       = ch_versions                 // channel: [ path(versions.yml) ]
//     //     if (!params.alignment_based)

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
