process PARSNP {
    // tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/parsnp:2.0.5--hdcf5f25_0' :
        'biocontainers/parsnp:2.0.5--hdcf5f25_0' }"

    input:
    path(contigs)
    path fasta

    output:
    tuple val(meta), path('*.xmfa')                              , emit: core_genome_alignment
    tuple val(meta), path('*.ggr')                               , emit: gingr_file
    tuple val(meta), path( "parsnp_output/parsnp.snps.mblocks" ) , emit: mblocks
    tuple val(meta), path( "parsnp_output/parsnp.tree" )         , emit: tree
    tuple val(meta), path( "parsnp_output/versions.yml" )        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parsnp \\
        -r $fasta \\
        -d $contigs \\
        -o ./parsnp_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version | cut -d ' ' -f 2 | sed 's/v//')
    END_VERSIONS
    """
}