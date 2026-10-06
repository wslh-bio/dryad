//
// Alignment_free subworkflow
// 

//
// Loading alignment free modules
//
include { MASHTREE                   } from '../../modules/local/mashtree'
include { ALIGNMENT_FREE_RESULTS     } from '../../modules/local/alignment_free_results'

//
// Creating alignment free workflow
//

workflow ALIGNMENT_FREE {

    take:
    reads               // channel: [ val(meta), [ reads ] ]
    cpus                // how many cpus mashtree should use
    quast_tsv           // will use quast summary output in final summary
    
    main:
    ch_versions = Channel.empty()       // Creating empty version channel to get versions.yml

    //
    // MASHTREE input: tuple val(meta), path(seqs)
    //
    MASHTREE (
        reads,
        cpus
    )

    ch_tree = MASHTREE.out.tree
    ch_versions = ch_versions.mix(MASHTREE.out.versions)

    if (!params.alignment_based) {
        ALIGNMENT_FREE_RESULTS (
        quast_tsv,
        params.run_name
        )

        ch_summary = ALIGNMENT_BASED_RESULTS.out.summary
    }
    else {
        ch_summary = channel.empty()
    }

    emit:
    tree           = ch_tree
    summary        = ch_summary
    versions       = ch_versions

}