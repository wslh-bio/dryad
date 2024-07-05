process IQTREE {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "staphb/iqtree2:2.3.4"

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
    iqtree2 \\
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
