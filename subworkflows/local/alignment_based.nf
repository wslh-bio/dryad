// Alignment_based subworkflow

include { REMOVE_REFERENCE } from '../../modules/local/remove_reference'
include { SAMPLE_COUNT     } from '../../modules/local/sample_count'
include { PARSNP           } from '../../modules/local/parsnp'
include { IQTREE           } from '../../modules/local/iqtree'
include { SNPDISTS         } from '../../modules/local/snpdists'
include { COMPARE_IO       } from '../../modules/local/compare_io'

workflow ALIGNMENT_BASED {

    take:
    reads               // channel: [ path[ reads ] ]
    fasta               // channel: /path/to/genome.fasta
    outdir              // output directory
    partition           // tells parsnp if it's important to partition
    add_reference       // tells parsnp if it needs to remove the reference
    samplesheet         // valid samplesheet to compare output to

    main:
    ch_versions = Channel.empty()       // Creating empty version channel to get versions.yml

//
// PARSNP
//
    PARSNP (
        reads,
        fasta,
        partition
        )
    ch_versions = ch_versions.mix(PARSNP.out.versions) 

//
// COMPARE_IO
//

    COMPARE_IO (
        samplesheet,
        PARSNP.out.tree
    )

//
// Remove reference
//
    if (!add_reference) {

        REMOVE_REFERENCE (
            PARSNP.out.mblocks
        )
        .set{ ch_for_mblocks }

        //
        // SAMPLE COUNT
        //
        SAMPLE_COUNT (
            ch_for_mblocks
        )

        //
        // IQTREE
        //
        IQTREE (
            ch_for_mblocks,
            SAMPLE_COUNT.out.count
        )
        ch_versions = ch_versions.mix(IQTREE.out.versions)

        //
        // SNPDISTS
        //
        SNPDISTS (
            ch_for_mblocks
        )
        ch_versions = ch_versions.mix(SNPDISTS.out.versions)
}

//
// Keep reference
//
    if (add_reference) {

    //
    // SAMPLE COUNT
    //
        SAMPLE_COUNT (
            PARSNP.out.mblocks
        )

    //
    // IQTREE
    //
        IQTREE (
            PARSNP.out.mblocks,
            SAMPLE_COUNT.out.count
        )
        ch_versions = ch_versions.mix(IQTREE.out.versions)

    //
    // SNPDISTS
    //
        SNPDISTS (
            PARSNP.out.mblocks
        )
        ch_versions = ch_versions.mix(SNPDISTS.out.versions)
    }

    emit:
    phylogeny    =      IQTREE.out.phylogeny
    tsv          =      SNPDISTS.out.tsv
    excluded     =      COMPARE_IO.out.excluded
    versions     =      ch_versions
}