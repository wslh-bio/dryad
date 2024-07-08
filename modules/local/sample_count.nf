process SAMPLE_COUNT {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    path fastas

    output:
    path "count.txt"   , emit: count

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    sample_count.py \\
    $fastas
    """
}