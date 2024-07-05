process MASHTREE {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "staphb/mashtree:1.4.6"

    input:
    path(seqs)
    val task_cpus

    output:
    path("*.dnd")           , emit: tree
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mashtree_bootstrap.pl --reps 100 --numcpus $task_cpus $seqs -- --min-depth 0 > mashtree.bootstrap.dnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
