process IQTREE {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "staphb/iqtree2:2.3.4"

    input:
    path(mblocks)
    path sample_count

    output:
    path("*.treefile")      , emit: phylogeny
    path("*.ufboot")        , emit: bootstrap, optional: true
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    number=\$(cat $sample_count)
    if [[ "\$number" -ge 4 ]]; then
        iqtree2 \\
                -s $mblocks \\
                -nt AUTO \\
                -m ${params.cg_tree_model} \\
                -bb 1000
    else
        iqtree2 \\
                -s $mblocks \\
                -nt AUTO \\
                -m ${params.cg_tree_model}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
