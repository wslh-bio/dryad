process SNPDISTS {

    label 'process_low'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "staphb/snp-dists:0.8.2"

    input:
    path(mblocks)

    output:
    path("snp_dists_matrix.tsv")            , emit: tsv
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    snp-dists \\
            -b $mblocks > snp_dists_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpdists: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
    """
}
