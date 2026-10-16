process PARSNP {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/wslh-bioinformatics/parsnp:2.1.6"

    input:
    path reads
    path fasta
    val partition
    val remove_recombination

    output:
    path( "parsnp_output/parsnp.xmfa"           )   , emit: core_genome_alignment
    path( "parsnp_output/parsnp.ggr"            )   , emit: gingr_file
    path( "parsnp_output/parsnp.snps.mblocks"   )   , emit: mblocks
    path( "parsnp_output/parsnp.tree"           )   , emit: tree, optional: true
    path( "parsnp_output/log/parsnpAligner.log" )   , emit: log
    path( "parsnp_output/*"                     )   , emit: all_output
    path( "versions.yml"                        )   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fasta_input = fasta.name != 'NO_FILE' ? "$fasta" : '!'
    def recombination = remove_recombination ? "-x" : ""
    """
    parsnp -r $fasta_input \\
            -d $reads \\
            -o ./parsnp_output \\
            $recombination \\
            $partition

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version | cut -d ' ' -f 2 | sed 's/v//')
    END_VERSIONS
    """

    stub:
    """
    mkdir parsnp_output
    mkdir parsnp_output/log
    touch "parsnp_output/parsnp.xmfa"
    touch "parsnp_output/parsnp.ggr"
    touch "parsnp_output/parsnp.snps.mblocks"
    touch "parsnp_output/parsnp.tree"
    touch "parsnp_output/log/parsnpAligner.log"
    touch "parsnp_output/*"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version | cut -d ' ' -f 2 | sed 's/v//')
    END_VERSIONS
    """
}