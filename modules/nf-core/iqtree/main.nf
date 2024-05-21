process IQTREE {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.3.0--h21ec9f0_0' :
        'biocontainers/iqtree:2.3.0--h21ec9f0_0' }"

    input:
    path(mblocks)

    output:
    path("*.treefile")      , emit: phylogeny
    path("*.ufboot")        , emit: bootstrap, optional: true
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    iqtree \\
            -s $mblocks \\
            -nt AUTO \\
            -m ${params.cg_tree_model} \\
            -bb 1000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
