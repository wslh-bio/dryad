process PARSE_PARSNP_ALIGNER_LOG {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path ("parsnp_output/log/parsnpAligner.log")
    val keep_ref

    output:
    path("aligner_log.tsv") , emit: aligner_log

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_parsnp_aligner_log.py parsnp_output/log/parsnpAligner.log $keep_ref
    """
}