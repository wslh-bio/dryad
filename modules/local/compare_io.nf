process COMPARE_IO {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path samplesheet
    path aligner_log

    output:
    path("excluded_samples_from_parsnp.txt"), emit: excluded

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin
    """
    compare_io.py \\
            $samplesheet \\ 
            $aligner_log
    """

}