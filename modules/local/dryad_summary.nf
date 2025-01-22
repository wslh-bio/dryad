process DRYAD_SUMMARY {

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
        path quast
        path aligner_log
        path excluded_samples

    output:
        path("*.csv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin
    """
    dryad_summary.py $aligner_log $quast $excluded_samples ${workflow.manifest.version}
    """

}