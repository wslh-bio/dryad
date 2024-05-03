//
// Alignment_free subworkflow
// 

//
// Loading alignment free modules
//
include { MASHTREE } from '../../modules/nf-core/mashtree'

//
// Creating alignment free workflow
//

workflow ALIGNMENT_FREE {

    take:
    reads       // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()       // Creating empty version channel to get versions.yml

    //
    // MASHTREE input: tuple val(meta), path(seqs)
    //
    MASHTREE (
        ch_samplesheet
    )

    //
    // Custom python script
    //

    emit:
    tree        = MASHTREE.out.tree
    matrix      = MASHTREE.out.matrix
    versions    = ch_versions

}