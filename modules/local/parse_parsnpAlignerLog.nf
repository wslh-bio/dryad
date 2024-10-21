PARSE_PARSNP_ALIGNER_LOG {

    input:
    path("parsnpAligner.log")

    output:
    path("")

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    """
}