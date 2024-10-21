PARSE_PARSNP_ALIGNER_LOG {

    input:
    path("parsnpAligner.log")

    output:
    path log
    val keep_ref

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_parsnpAlignerLog.py $log $keep_ref
    """
}