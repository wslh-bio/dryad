process REMOVE_REFERENCE {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    path fastas
    path reference_fasta


    output:
    path( "parsnp_output/cleaned_parsnp.snps.mblocks"   )   , emit: mblocks

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    remove_reference.py \\
    $fastas \\
    $reference_fasta
    """
}