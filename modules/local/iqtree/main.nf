process IQTREE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "staphb/iqtree2:2.3.4"

    input:
    path(mblocks)

    output:
    path("*.treefile")      , emit: phylogeny, optional: true
    path("*.ufboot")        , emit: bootstrap, optional: true
    path "versions.yml"     , emit: versions
    path "note.txt"         , emit: note, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sample_count=\$( grep -c '>' $mblocks)

    if [[ \$sample_count -ge 4 ]]; then
        iqtree2 \\
                -s $mblocks \\
                -nt AUTO \\
                -m ${params.cg_tree_model} \\
                -bb 1000
    elif [[ \$sample_count -eq 3  ]]; then
        iqtree2 \\
                -s $mblocks \\
                -nt AUTO \\
                -m ${params.cg_tree_model}
    else
        echo "IQ-TREE was not run, because at least three samples were not present in the alignment." > note.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree2: \$(echo \$(iqtree2 -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch "cleaned_parsnp.snps.mblocks.treefile"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree2: \$(echo \$(iqtree2 -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
