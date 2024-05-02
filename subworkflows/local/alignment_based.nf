// Alignment_based subworkflow
// Parsnp>IQ-TREE&snp-dists>python
include { PARSNP   }      from '../modules/nf-core/parsnp'
include { IQTREE   }      from '../modules/nf-core/iqtree'
include { SNPDISTS }      from '../modules/nf-core/snpdists'

workflow ALIGNMENT_BASED {

    take:
    reads       // channel: [ val(meta), [ reads ] ]
    fasta       // channel: /path/to/genome.fasta

    main:
    ch_versions = Channel.empty()

//
// PARSNP
//
    PARSNP (
        reads,
        fasta
        )
    ch_versions = ch_versions.mix(PARSNP.out.versions.first()) 

//
// Creating channel for phylogeny to be fed into IQTREE and 
//
    PARSNP
        .out
        .phylogeny
        .set { ch_phylogeny }

    PARSNP
        .out
        .mblocks
        .set { ch_mblocks }

    IQTREE (
    ch_phylogeny
    )
    ch_versions = ch_versions.mix(IQTREE.out.versions.first())

    SNPDISTS (
        ch_mblocks
    )
    ch_versions = ch_versions.mix(SNPDISTS.out.versions.first())

    emit:
    tuple val(meta), path("*.treefile"),    emit: phylogeny
    tuple val(meta), path("*.tsv")     ,    emit: tsv
    path "versions.yml"                ,    emit: versions
}