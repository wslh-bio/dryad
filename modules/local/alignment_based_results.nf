process ALIGNMENT_BASED_RESULTS {

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
        path quast
        path aligner_log
        path excluded_samples
        val run_name

    output:
        path("*.csv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin
    def cleaned_runname=run_name.toString().replaceAll(' ', '_')
    """
    summarize_alignment_based.py \\
        --quast_results $quast \\
        --run_name $cleaned_runname \\
        --dryad_version ${workflow.manifest.version} \\
        --aligner_log $aligner_log \\
        --excluded_samples $excluded_samples
    """

}