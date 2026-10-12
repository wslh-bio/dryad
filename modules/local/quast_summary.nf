process QUAST_SUMMARY {
    label 'process_single'
    
    container "quay.io/wslh-bioinformatics/pandas:1.5.0"

    input:
    path("data/*")

    output:
    path("quast_results.tsv")    , emit: quast_tsv
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quast_summary.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$( python3 --version 2>&1 ) | sed 's/^.*Python //' )
    END_VERSIONS
    """
}