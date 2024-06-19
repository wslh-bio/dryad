process SNPDISTS {

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snp-dists:0.8.2--h5bf99c6_0' :
        'biocontainers/snp-dists:0.8.2--h5bf99c6_0' }"

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
