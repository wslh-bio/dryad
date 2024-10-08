process PARSNP {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "staphb/parsnp:2.0.5"

    input:
    path reads
    path fasta
    val partition

    output:
    path( "parsnp_output/parsnp.xmfa"           )   , emit: core_genome_alignment
    path( "parsnp_output/parsnp.ggr"            )   , emit: gingr_file
    path( "parsnp_output/parsnp.snps.mblocks"   )   , emit: mblocks
    path( "parsnp_output/parsnp.tree"           )   , emit: tree
    path( "versions.yml"                        )   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parsnp -r $fasta \\
           -d $reads \\
           -o ./parsnp_output \\
           $partition

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version | cut -d ' ' -f 2 | sed 's/v//')
    END_VERSIONS
    """
}