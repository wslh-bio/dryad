process COMPARE_IO {
    label 'process_single'
    
    container "quay.io/wslh-bioinformatics/pandas:1.5.0"
    
    input:
    path samplesheet
    path aligner_log
    val fasta

    output:
    path("sample_exclusion_status.csv")    , emit: excluded
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin
    def fasta_input = fasta.name != 'NO_FILE' ? "$fasta" : 'random'
    """
    compare_io.py \\
    --input_file $samplesheet \\
    --output_file $aligner_log \\
    --reference_genome $fasta_input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$( python3 --version 2>&1 ) | sed 's/^.*Python //' )
    END_VERSIONS
    """

    stub:
    """
    touch "sample_exclusion_status.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$( python3 --version 2>&1 ) | sed 's/^.*Python //' )
    END_VERSIONS
    """

}