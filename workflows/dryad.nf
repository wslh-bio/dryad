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
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? channel.fromPath( params.multiqc_config, checkIfExists: true ) : channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo, checkIfExists: true ) : channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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
include { INPUT_CHECK                    } from '../subworkflows/local/input_check'
include { COUNT_FASTA                    } from '../modules/local/count_fasta'
include { REJECTED_SAMPLES               } from '../modules/local/rejected_samples'
include { QUAST                          } from '../modules/local/quast'
include { QUAST_SUMMARY                  } from '../modules/local/quast_summary'
include { ALIGNMENT_BASED                } from '../subworkflows/local/alignment_based'
include { ALIGNMENT_FREE                 } from '../subworkflows/local/alignment_free'
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                        } from '../modules/local/multiqc/multiqc'

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
    .set { ch_input }

    // Adding version information
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Run Module: count_fasta
    COUNT_FASTA(
        ch_input
    )

    ch_csv = COUNT_FASTA.out.csv
                .splitCsv(header: true)
                .join(ch_input)
                .map { meta, csv, file ->
                def count = csv.count as Integer
                tuple(meta, file, count)
                }

    ch_versions = ch_versions.mix(COUNT_FASTA.out.versions)

    // Pass/fail based on sequence (contig) count in fasta files
    ch_csv
        .branch{ meta, file, count ->
            pass: count > params.seqcount_cutoff
            fail: count <= params.seqcount_cutoff
        }
        .set{ ch_fasta }

    ch_fasta.pass
        .map { meta, file, _count -> 
            [meta, file]
            }
        .set{ ch_filtered }

    ch_fasta.fail
        .map { meta, file, _count ->
            meta.id
            }
        .set{ ch_failed }

    // Collect empty/rejected samples (if any)
    ch_failed
        .ifEmpty('NO_EMPTY_SAMPLES')
        .collectFile(
                name: 'empty_samples.csv',
                newLine: true
            )
        .set{ ch_rejected_file }

    REJECTED_SAMPLES (
        ch_rejected_file,
        "dryad",
        params.run_name
    )

    ch_versions = ch_versions.mix(REJECTED_SAMPLES.out.versions)

    //
    // QC check for runs if skip quast
    //
    if (!params.skip_quast) {
        QUAST ( ch_filtered )
        ch_quast_reports = QUAST.out.transposed_report
        ch_versions = ch_versions.mix(QUAST.out.versions)

        QUAST_SUMMARY (
            ch_quast_reports.collect()
            )
        ch_quast_summary = QUAST_SUMMARY.out.quast_tsv
        ch_versions = ch_versions.mix(QUAST_SUMMARY.out.versions)

    } else {
        quast_file = file("$baseDir/assets/empty.txt",checkIfExists:true)
    }

    //
    // Re-mapping channel to intake paths
    //
    ch_filtered
        .map { sample, fasta ->
            fasta
        } // Produces queue channel of just fasta file paths in a list
        .collect()
        .set { ch_fasta_paths }

    //
    // SUBWORKFLOW: Alignment Free
    //
    if (params.alignment_free) {
        if (!params.skip_quast) {
            ALIGNMENT_FREE (
                ch_fasta_paths,
                params.task.cpus,
                ch_quast_summary
                )
        }
        if (params.skip_quast) {
            ALIGNMENT_FREE (
                ch_fasta_paths,
                params.task.cpus,
                quast_file
                )
        }
        ch_alignment_free_versions = ALIGNMENT_FREE.out.versions
        ch_versions = ch_versions.mix(ch_alignment_free_versions)
    }

    if (params.fasta == "random") {
        ch_fasta = file(params.random_file, checkIfExists:true)
    }
    else if (!params.fasta ) {
        ch_fasta = Channel.empty()
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
                ch_fasta_paths,
                ch_fasta,
                params.outdir,
                params.parsnp_partition,
                params.add_reference,
                params.remove_recombination,
                INPUT_CHECK.out.csv,
                QUAST_SUMMARY.out.quast_tsv
                )
            ch_alignment_based_versions = ALIGNMENT_BASED.out.versions
            ch_versions = ch_versions.mix(ch_alignment_based_versions)
        }
        if (params.skip_quast) {
            ALIGNMENT_BASED (
                ch_fasta_paths,
                ch_fasta,
                params.outdir,
                params.parsnp_partition,
                params.add_reference,
                params.remove_recombination,
                INPUT_CHECK.out.csv,
                quast_file
                )
            ch_alignment_based_versions = ALIGNMENT_BASED.out.versions
            //ch_versions = ch_versions.mix(ch_alignment_based_versions)
        }        
    }
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    //MODULE: MultiQC
    //
    
    workflow_summary    = WorkflowDryad.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = channel.value(workflow_summary)

    methods_description    = WorkflowDryad.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = channel.value(methods_description)

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.result.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
}