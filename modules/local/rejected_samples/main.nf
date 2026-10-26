process REJECTED_SAMPLES {

    label 'process_low'
    container "quay.io/nf-core/ubuntu:20.04"

    input:
    path csv_file
    val workflow_name
    val run_name

    output:
    path "*_rejected_samples.csv", emit: rejected_samples_report
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash
    cp ${csv_file} ${run_name}_${workflow_name}_rejected_samples.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$( bash -c 'echo \$BASH_VERSION' )
    END_VERSIONS
    """

    stub:
    """
    touch "${workflow_name}_empty_samples.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$( bash -c 'echo \$BASH_VERSION' )
    END_VERSIONS
    """
}
