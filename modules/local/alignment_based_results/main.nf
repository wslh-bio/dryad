process CREATE_REPORT {

    label 'process_single'
    container "quay.io/wslh-bioinformatics/pandas:1.5.0"
    
    input:
    path quast
    path aligner_log
    path excluded_samples
    val run_name

    output:
    path("*.csv"), emit: summary
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_report: \$( echo \$( python3 --version 2>&1 ) | sed 's/^.*Python //' )
    END_VERSIONS
    """
}