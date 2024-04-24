//
// Alignment_free subworkflow
// Mashtree

workflow ALIGNMENT_FREE {

    take:
        ch_input_reads      // channel: [ val(meta), path(input)]

    main:
        ch_versions = Channel.empty()
        ch_input    = ch_input_reads

    emit:
        tuple val(meta), path("*.dnd"), emit: tree
        tuple val(meta), path("*.tsv"), emit: matrix
        path "versions.yml"           , emit: versions

    //
    // MASHTREE input: tuple val(meta), path(seqs)
    //
    MASHTREE (
        ch_input
    )

    //
    // Custom python script
    //

    
}