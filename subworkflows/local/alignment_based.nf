// Alignment_based subworkflow

include { REMOVE_REFERENCE } from '../../modules/local/remove_reference'
include { SAMPLE_COUNT               } from '../../modules/local/sample_count'
include { PARSNP                     } from '../../modules/local/parsnp'
include { IQTREE                     } from '../../modules/local/iqtree'
include { SNPDISTS                   } from '../../modules/local/snpdists'
include { PARSE_PARSNP_ALIGNER_LOG   } from '../../modules/local/parse_parsnp_aligner_log'
include { COMPARE_IO                 } from '../../modules/local/compare_io'
include { DRYAD_SUMMARY              } from '../../modules/local/dryad_summary'

workflow ALIGNMENT_BASED {

    take:
    reads               // channel: [ path[ reads ] ]
    fasta               // channel: /path/to/genome.fasta
    outdir              // output directory
    partition           // tells parsnp if it's important to partition
    add_reference       // tells parsnp if it needs to remove the reference
    samplesheet         // valid samplesheet to compare output to
    quast_tsv           // will use quast summary output in final summary

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
        // PARSER
        //
        PARSE_PARSNP_ALIGNER_LOG (
            PARSNP.out.log,
            add_reference
        )

        //
        // COMPARE_IO
        //
        COMPARE_IO (
            samplesheet,
            PARSE_PARSNP_ALIGNER_LOG.out.aligner_log
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
        
    //
    // Final Summary
    //
    DRYAD_SUMMARY (
        quast_tsv,
        PARSE_PARSNP_ALIGNER_LOG.out.aligner_log,
        COMPARE_IO.out.excluded
        )
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
    // PARSER
    //
    PARSE_PARSNP_ALIGNER_LOG (
        PARSNP.out.log,
        add_reference
    )

    //
    // COMPARE_IO
    //
    COMPARE_IO (
        samplesheet,
        PARSE_PARSNP_ALIGNER_LOG.out.aligner_log
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
    
    //
    // Final Summary
    //
    DRYAD_SUMMARY (
        quast_tsv,
        PARSE_PARSNP_ALIGNER_LOG.out.aligner_log,
        COMPARE_IO.out.excluded
        )
    }

    emit:
    phylogeny    =      IQTREE.out.phylogeny
    tsv          =      SNPDISTS.out.tsv
    aligner_log  =      PARSE_PARSNP_ALIGNER_LOG.out.aligner_log
    excluded     =      COMPARE_IO.out.excluded
    summary      =      DRYAD_SUMMARY.out.summary
    versions     =      ch_versions
}