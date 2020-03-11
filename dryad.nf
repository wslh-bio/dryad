#!/usr/bin/env nextflow

//Description: Workflow for building various trees from raw illumina reads
//Author: Kelsey Florek
//eMail: kelsey.florek@slh.wisc.edu

//starting parameters
params.reads = ""
params.outdir = "dryad_results"
params.cg = false
params.snp = false
params.snp_reference = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.fastq.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
    .set { raw_reads }

Channel
    .fromPath(params.snp_reference)
    .set { snp_reference }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(oldName), file(reads) from raw_reads

  output:
  tuple name, file("${name}_{R1,R2}.fastq.gz") into read_files_fastqc, read_files_trimming

  script:
  name = oldName.split("_")[0]
  """
  mv ${reads[0]} ${name}_R1.fastq.gz
  mv ${reads[1]} ${name}_R2.fastq.gz
  """
}

//Step1a: FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from read_files_fastqc

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"
  publishDir "${params.outdir}/trimmed", mode: 'copy'

  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}*{_1,_2}.fastq.gz") into trimmed_reads
  file("${name}.trim.stats.txt") into trimmed_reads_stats

  script:
  """
  cpus=`grep -c ^processor /proc/cpuinfo`
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads \$cpus ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} > ${name}.trim.stats.txt
  mv ${name}*1P.fastq.gz ${name}_1.fastq.gz
  mv ${name}*2P.fastq.gz ${name}_2.fastq.gz
  """
}
//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"
  publishDir "${params.outdir}/cleanedreads", mode: 'copy',pattern:"*.stats.txt"

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads_cg,cleaned_reads_snp
  file("${name}.{phix,adapters}.stats.txt") into read_cleanning_stats

  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}

//SNP Step1: Run CFSAN-SNP Pipeline
process cfsan {
  publishDir "${params.outdir}/cfsan-snp", mode: 'copy'

  input:
  val(items) from cleaned_reads_snp.toList().collect()
  file(reference) from snp_reference

  output:
  file("snp_distance_matrix.tsv")
  file("snpma.fasta") into snp_alignment

  when:
  params.snp == true

  script:
  //create input reads dir in work dir
  new File("input_reads").mkdir()
  items.each{
    //iterate over the files and get the name and move reads to subfolders with name
    String name = it.get(0)
    new File("input_reads/$name").mkdir()
    def read1 = it.get(1).get(0)
    def read2 = it.get(1).get(1)
    def read1_dest = new File("input_reads/$name/",read1.name)
    read1.withInputStream{stream-> read1_dest << stream }
    def read2_dest = new File("input_reads/$name/",read2.name)
    read2.withInputStream{stream-> read2_dest << stream }
  }
  """
  cfsan_snp_pipeline run ${reference} -o . -s input_reads
  """

}

//SNP Step2: Run IQTREE on snp alignment
process snp_tree {
  publishDir "${params.outdir}/snp_tree", mode: 'copy'

  input:
  file(snp_fasta) from snp_alignment

  output:
  file("snp.tree")

  script:
  numGenomes = 0
  snp_fasta.eachLine {
    line -> if(line.startsWith('>')){numGenomes += 1}
  }
  if(numGenomes > 3){
    """
    iqtree -nt AUTO -s core_gene_alignment.aln -m ${params.cg_tree_model} -bb 1000
    mv snpma.fasta.contree snp.tree
    """
  }
  else {
    """
    iqtree -nt AUTO -s core_gene_alignment.aln -m ${params.cg_tree_model}
    mv snpma.fasta.treefile snp.tree
    """
  }
}

//CG Step1: Assemble trimmed reads with Shovill
process shovill {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/assembled", mode: 'copy'

  input:
  set val(name), file(reads) from cleaned_reads_cg

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_annotation

  shell:
  '''
  ram=`awk '/MemAvailable/ { printf "%.0f \\n", $2/1024/1024 }' /proc/meminfo`
  shovill --cpus 0 --ram $ram  --outdir . --R1 !{reads[0]} --R2 !{reads[1]} --force
  mv contigs.fa !{name}.contigs.fa
  '''
}

//CG Step2a: Assembly Quality Report
process quast {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  file(assemblies) from assembled_genomes_quality.collect()

  output:
  file "assembly.report.txt" into quast_report

  script:
  """
  quast.py ${assemblies} -o .
  mv report.txt assembly.report.txt
  """
}


//CG Step2b: Annotate with prokka
process prokka {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/annotated",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_annotation

  output:
  file "${name}.gff" into annotated_genomes

  when:
  params.cg == true

  script:
  """
  prokka --cpu 0 --force --compliant --prefix ${name} --mincontiglen 500 --outdir . ${assembly}
  """
}

//CG Step3: Align with Roary
process roary {
  publishDir "${params.outdir}/core_alignment",mode:'copy'
  numGenomes = 0
  input:
  file(genomes) from annotated_genomes.collect()

  output:
  file("core_gene_alignment.aln") into core_aligned_genomes
  file "core_genome_statistics.txt" into core_aligned_stats

  script:
  if(params.roary_mafft == true){
    mafft="-n"
  }else{mafft=""}
  """
  cpus=`grep -c ^processor /proc/cpuinfo`
  roary -e ${mafft} -p \$cpus ${genomes}
  mv summary_statistics.txt core_genome_statistics.txt
  """
}

//CG Step4: IQTree for core-genome
process cg_tree {
  publishDir "${params.outdir}/core_genome_tree",mode:'copy'

  input:
  file(alignedGenomes) from core_aligned_genomes

  output:
  file("core_genome.tree")


  script:
  numGenomes = 0
  alignedGenomes.eachLine {
    line -> if(line.startsWith('>')){numGenomes += 1}
  }
  if(numGenomes > 3){
    """
    iqtree -nt AUTO -s core_gene_alignment.aln -m ${params.cg_tree_model} -bb 1000
    mv core_gene_alignment.aln.contree core_genome.tree
    """
  }
  else {
    """
    iqtree -nt AUTO -s core_gene_alignment.aln -m ${params.cg_tree_model}
    mv core_gene_alignment.aln.treefile core_genome.tree
    """
  }
}

//Collect Results
process multiqc {
  tag "$prefix"
  publishDir "${params.outdir}/MultiQC", mode: 'copy'
  echo true

  input:
  //file multiqc_config
  file (fastqc:'fastqc/*') from fastqc_results.collect()
  file ('*.trim.stats.txt') from trimmed_reads_stats.collect()
  file ('assembly.report.txt') from quast_report
  file ('core_genome_statistics.txt') from core_aligned_stats
  file ('*.stats.txt') from read_cleanning_stats.collect()


  output:
  file "*multiqc_report.html" into multiqc_report
  file "*_data"

  script:
  prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
  """
  multiqc . 2>&1
  """
}
