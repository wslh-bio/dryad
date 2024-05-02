//
// Alignment_free subworkflow
// 

//
// Loading alignment free modules
//
include { MASHTREE } from '../modules/nf-core/mashtree'

//
// Creating alignment free workflow
//

workflow ALIGNMENT_FREE {

    ch_versions = Channel.empty()

    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    INPUT_CHECK.out.reads
        .set{ ch_input_reads}

    //
    // MODULE: MashTree
    //
    MASHTREE (
        ch_input 
    )
    .map {
            meta, fasta ->
                single : fasta.size() == 1
                    return [ meta, fasta.flatten() ]
                multiple: fasta.size() > 1
                    return [ meta, fasta.flatten() ]
    }

    ch_versions = ch_versions.mix(MASHTREE.out.versions.first().ifEmpty(null))

    // take:
    //     ch_samplesheet      // channel: [ val(meta), path(input)]

    // main:
    //     ch_versions = Channel.empty()
    //     ch_input    = ch_input_reads

    // emit:
    //     tuple val(meta), path("*.dnd"), emit: tree
    //     tuple val(meta), path("*.tsv"), emit: matrix
    //     path "versions.yml"           , emit: versions

    // //
    // // MASHTREE input: tuple val(meta), path(seqs)
    // //
    // MASHTREE (
    //     ch_samplesheet
    // )

    //
    // Custom python script
    //

}