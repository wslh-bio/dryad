process COMPARE_IO {

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path samplesheet
    path parsnp_tree

    output:
    path "excluded_samples_from_parsnp.txt"   , emit: excluded

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin
    """
    compare_io.py \\
            $samplesheet \\ 
            $parsnp_tree
    """

}