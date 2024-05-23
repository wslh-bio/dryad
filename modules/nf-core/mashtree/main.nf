process MASHTREE {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mashtree:1.4.6--pl5321h031d066_0' :
        'biocontainers/mashtree:1.4.6--pl5321h031d066_0' }"

    input:
    path(seqs)

    output:
    path("*.dnd")           , emit: tree
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mashtree_bootstrap.pl --reps 100 --numcpus ${task.cpus} $seqs -- --min-depth 0 > mashtree.bootstrap.dnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
