// Alignment_based subworkflow
// Parsnp>IQ-TREE&snp-dists>python
if (params.reference) { ch_reference = file(params.reference) } else { exit 1, 'Reference fasta file not specified!' }

workflow ALIGNMENT_BASED {

    take:
        ch_samplesheet

    main:
        ch_versions = Channel.empty()
        // ch_input    = ch_input_reads
    
    emit:
        tuple val(meta), path("*.dnd"), emit: tree
        tuple val(meta), path("*.tsv"), emit: matrix
        path "versions.yml"           , emit: versions

    MASHTREE (
        ch_input
    )
}