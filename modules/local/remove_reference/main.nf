process REMOVE_REFERENCE {

    container "quay.io/biocontainers/biopython:1.81"

    input:
    path fastas

    output:
    path "*.mblocks"        , emit: mblocks
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    remove_reference.py \\
    $fastas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$( python3 --version 2>&1 ) | sed 's/^.*Python //' )
    END_VERSIONS
    """

    stub:
    """
    touch "cleaned_parsnp.snps.mblocks"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$( python3 --version 2>&1 ) | sed 's/^.*Python //' )
    END_VERSIONS
    """
}