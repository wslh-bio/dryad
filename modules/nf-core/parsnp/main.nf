process PARSNP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/parsnp:2.0.5--hdcf5f25_0' }"

    input:
    tuple val(meta), path(reads)
    path fasta
    path outdir

    output:
    path "*.xmfa"               , emit: core-genome-alignment
    path "*.ggr"                , emit: gingr-file
    path "parsnp.snps.mblocks"  , emit: mblocks
    path "parsnp.tree"          , emit: phylogeny
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parsnp \\
        -r ${fasta} \\
        -d ${reads} \\
        -o ${outdir}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            parsnp: \$(parsnp -V | cut -d ' ' -f2 | sed 's/v//')
        END_VERSIONS
    """
}