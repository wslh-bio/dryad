process QUAST_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("data/*")

    output:
    path("quast_results.tsv"), emit: quast_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quast_summary.py
    """
}