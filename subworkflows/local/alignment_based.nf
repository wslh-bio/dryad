// Alignment_based subworkflow

include { REMOVE_REFERENCE           } from '../../modules/local/remove_reference'
include { PARSNP                     } from '../../modules/local/parsnp'
include { IQTREE                     } from '../../modules/local/iqtree'
include { SNPDISTS                   } from '../../modules/local/snpdists'
include { PARSE_ALIGNER_LOG          } from '../../modules/local/parse_parsnp_aligner_log'
include { COMPARE_IO                 } from '../../modules/local/compare_io'
include { ALIGNMENT_BASED_RESULTS    } from '../../modules/local/alignment_based_results'

workflow ALIGNMENT_BASED {

    take:
    reads               // channel: [ path[ reads ] ]
    fasta               // channel: /path/to/genome.fasta
    outdir              // output directory
    partition           // tells parsnp if it's important to partition
    add_reference       // tells parsnp if it needs to remove the reference
    recombination       // tells parsnp to remove_recombination
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
        partition,
        recombination
        )
    
    ch_mblocks = PARSNP.out.mblocks
    ch_parsnp_log = PARSNP.out.log
    ch_versions = ch_versions.mix(PARSNP.out.versions) 

//
// Remove reference
//
    if (!add_reference) {
        REMOVE_REFERENCE (
            ch_mblocks
        )
        .set{ ch_rmref_mblocks }

        //
        // PARSER
        //
        PARSE_ALIGNER_LOG (
            ch_parsnp_log,
            add_reference
        )

        ch_parsed = PARSE_ALIGNER_LOG.out.aligner_log

        //
        // COMPARE_IO
        //
        COMPARE_IO (
            samplesheet,
            ch_parsed
        )

        ch_excluded_samples = COMPARE_IO.out.excluded
        
        //
        // IQTREE
        //
        IQTREE (
            ch_rmref_mblocks
        )

        ch_tree = IQTREE.out.phylogeny
        ch_versions = ch_versions.mix(IQTREE.out.versions)

        //
        // SNPDISTS
        //
        SNPDISTS (
            ch_rmref_mblocks
        )

        ch_snp = SNPDISTS.out.tsv
        ch_versions = ch_versions.mix(SNPDISTS.out.versions)

        //
        // Final Summary
        //
        ALIGNMENT_BASED_RESULTS (
            quast_tsv,
            ch_parsed,
            ch_excluded_samples,
            params.run_name
            )
            
        ch_summary = ALIGNMENT_BASED_RESULTS.out.summary
    }

//
// Keep reference
//
    if (add_reference) {

        //
        // PARSER
        //
        PARSE_ALIGNER_LOG (
            ch_parsnp_log,
            add_reference
        )

        ch_parsed = PARSE_ALIGNER_LOG.out.aligner_log

        //
        // COMPARE_IO
        //
        COMPARE_IO (
            samplesheet,
            ch_parsed
        )

        ch_excluded_samples = COMPARE_IO.out.excluded

        //
        // IQTREE
        //
        IQTREE (
            PARSNP.out.mblocks
        )

        ch_tree = IQTREE.out.phylogeny
        ch_versions = ch_versions.mix(IQTREE.out.versions)

        //
        // SNPDISTS
        //
        SNPDISTS (
            PARSNP.out.mblocks
        )

        ch_snp = SNPDISTS.out.tsv
        ch_versions = ch_versions.mix(SNPDISTS.out.versions)

        //
        // Final Summary
        //
        ALIGNMENT_BASED_RESULTS (
            quast_tsv,
            ch_parsed,
            ch_excluded_samples,
            params.runname
            )

        ch_summary = ALIGNMENT_BASED_RESULTS.out.summary
    }

    emit:
    phylogeny    =      ch_tree
    tsv          =      ch_snp 
    aligner_log  =      ch_parsed
    excluded     =      ch_excluded_samples
    summary      =      ch_summary
    versions     =      ch_versions
}
