#!/usr/bin/env nextflow

//Description: Workflow for generating a genomic comparison between HAI/AR samples.
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

params.test = false

if(params.test){
  testIDS = ['SRR14311557','SRR14311556','SRR14311555','SRR14311554',
    'SRR14311553','SRR14311552','SRR14613509','SRR14874874']
  println "Running test analysis using the following samples:"
  println testIDS
  Channel
      .fromSRA(testIDS)
      .into { raw_reads; raw_reads_count }
  params.snp_reference = "$baseDir/assets/ASM211692v1.fasta"

} else{
  //setup channel to read in and pair the fastq files
  Channel
      .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
      .into { raw_reads; raw_reads_count }
  //check we have at least 3 samples
  Channel
      .from(raw_reads_count)
      .collect()
      .subscribe {
        int size = it.queue[0].size()
        if(size < 3){
          println "Dryad requires 3 or more samples."
          System.exit(1)
        }
      }
}

if (params.snp_reference) {
    Channel
        .fromPath(params.snp_reference)
        .set { snp_reference }
    Channel
        .from("$baseDir/snppipeline.conf")
        .set { snp_config }
}

//Step0: Preprocess reads - change names
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(outfiles) into read_files_fastqc, read_files_trimming, read_files_kraken

  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    outfiles = ["${name}_R1.fastq.gz","${name}_R2.fastq.gz"]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
    outfiles = reads
    """
    """
  }
}

//Step1: Trim reads and remove PhiX contamination
process clean_reads {
  tag "$name"
  publishDir "${params.outdir}/trimming", mode: 'copy',pattern:"*.trim.txt"

  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_shovill, cleaned_reads_fastqc
  file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_snp
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats
  tuple file("${name}.phix.stats.txt"),file("${name}.adapters.stats.txt"),file("${name}.trim.txt") into multiqc_clean_reads

  script:
  """
  bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.trimmed_1.fastq.gz out2=${name}.trimmed_2.fastq.gz qtrim=window,${params.windowsize} trimq=${params.qualitytrimscore} minlength=${params.minlength} tbo tbe &> ${name}.out
  repair.sh in1=${name}.trimmed_1.fastq.gz in2=${name}.trimmed_2.fastq.gz out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_clean_1.fastq.gz out2=${name}_clean_2.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${name}.out > ${name}.trim.txt
  """
}

combined_reads = read_files_fastqc.concat(cleaned_reads_fastqc)

//QC Step: FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from combined_reads

  output:
  file("*_fastqc.{zip,html}") into fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//QC Step: Summarize FastQC
process fastqc_summary {
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  file(fastqc) from fastqc_results.collect()

  output:
  file("fq_summary.txt") into fastqc_summary

  shell:
  """
  zips=`ls *.zip`
  for i in \$zips; do
      unzip -o \$i &>/dev/null;
  done
  fq_folders=\${zips}
  for folder in \$fq_folders; do
    folder=\${folder%.*}
    cat \$folder/summary.txt >> fq_summary.txt
    ls .
  done;
  sed -i 's/.fastq.gz//g' fq_summary.txt
  """
}

//CG Step1: Assemble trimmed reads with Shovill and map reads back to assembly
process shovill {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/assembled", mode: 'copy',pattern:"*.fa"
  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.sam"

  input:
  set val(name), file(reads) from cleaned_reads_shovill

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_annotation, assembled_genomes_quality
  tuple name, file("${name}.sam") into sam_files

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory} --outdir ./output --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv ./output/contigs.fa ${name}.contigs.fa
  bwa index ${name}.contigs.fa
  bwa mem ${name}.contigs.fa ${reads[0]} ${reads[1]} > ${name}.sam
  """
}

//QC Step: Index and sort bam file then calculate coverage
process samtools {
  tag "$name"

  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.bam"
  publishDir "${params.outdir}/coverage", mode: 'copy', pattern:"*_depth.tsv*"

  input:
  set val(name), file(sam) from sam_files

  output:
  file("${name}_depth.tsv") into cov_files
  file("${name}.stats.txt") into stats_multiqc

  shell:
  """
  samtools view -S -b ${name}.sam > ${name}.bam
  samtools sort ${name}.bam > ${name}.sorted.bam
  samtools index ${name}.sorted.bam
  samtools depth -a ${name}.sorted.bam > ${name}_depth.tsv
  samtools stats ${name}.sorted.bam > ${name}.stats.txt
  """
}

//QC Step: Calculate coverage stats
process coverage_stats {
  publishDir "${params.outdir}/coverage", mode: 'copy'

  input:
  file(cov) from cov_files.collect()

  output:
  file('coverage_stats.txt')

  script:
  """
  #!/usr/bin/env python3
  import glob
  import os
  from numpy import median
  from numpy import average

  results = []

  files = glob.glob("*_depth.tsv*")
  for file in files:
    nums = []
    sid = os.path.basename(file).split('_')[0]
    with open(file,'r') as inFile:
      for line in inFile:
        nums.append(int(line.strip().split()[2]))
      med = median(nums)
      avg = average(nums)
      results.append(f"{sid}\\t{med}\\t{avg}\\n")

  with open('coverage_stats.txt', 'w') as outFile:
    outFile.write("Sample\\tMedian Coverage\\tAverage Coverage\\n")
    for result in results:
      outFile.write(result)
  """
}
//QC Step: Run Quast on assemblies
process quast {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_multiqc

  script:
  """
  quast.py ${assembly} -o .
  mv report.txt ${name}.quast.tsv
  """
}

if (params.snp_reference) {
    //SNP Step1: Run CFSAN-SNP Pipeline
    process cfsan {
      publishDir "${params.outdir}", mode: 'copy'

      input:
      file(reads) from cleaned_reads_snp.collect()
      file(reference) from snp_reference
      file(config) from snp_config
      output:
      file("snp_distance_matrix.tsv") into snp_mat
      file("snpma.fasta") into snp_alignment

      script:
      """
      #!/usr/bin/env python
      import subprocess
      import glob
      import os,sys

      fwd_reads = glob.glob("*_clean_1.fastq.gz")
      fwd_reads.sort()
      rev_reads = glob.glob("*_clean_2.fastq.gz")
      rev_reads.sort()

      if len(fwd_reads) != len(rev_reads):
        sys.exit("Uneven number of forward and reverse reads.")

      os.mkdir("input_reads")

      c = 0
      while c < len(fwd_reads):
        name = os.path.basename(fwd_reads[c]).split('_clean_1')[0]
        path = os.path.join("input_reads",name)
        os.mkdir(path)
        os.rename(fwd_reads[c],os.path.join(path,fwd_reads[c]))
        os.rename(rev_reads[c],os.path.join(path,rev_reads[c]))
        c += 1

      command = "cfsan_snp_pipeline run ${reference} -o . -s input_reads"
      process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
      output, error = process.communicate()
      print output
      print error
      """
    }

    //SNP Step2: Run IQTREE on snp alignment
    process snp_tree {
      publishDir "${params.outdir}", mode: 'copy'

      input:
      file(snp_fasta) from snp_alignment

      output:
      file("snp.tree") optional true

      script:
        """
        numGenomes=`grep -o '>' snpma.fasta | wc -l`
        if [ \$numGenomes -gt 3 ]
        then
          iqtree -nt AUTO -s snpma.fasta -m ${params.cg_tree_model} -bb 1000
          mv snpma.fasta.contree snp.tree
        fi
        """
    }
}

//CG Step2: Annotate with prokka
//TODO: add genus and species
process prokka {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/annotated",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_annotation

  output:
  file("${name}.gff") into annotated_genomes
  file("${name}.prokka.stats.txt") into prokka_multiqc

  script:
  """
  prokka --cpu ${task.cpus} --force --compliant --prefix ${name} --mincontiglen 500 --outdir . ${assembly} > ${name}.log
  mv ${name}.txt ${name}.prokka.stats.txt
  """
}
//CG Step3: Align with Roary
process roary {
  publishDir "${params.outdir}",mode:'copy'

  numGenomes = 0
  input:
  file(genomes) from annotated_genomes.collect()

  output:
  file("core_gene_alignment.aln") into core_aligned_genomes
  file("core_genome_statistics.txt") into core_aligned_stats

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
  publishDir "${params.outdir}",mode:'copy'

  input:
  file(alignedGenomes) from core_aligned_genomes

  output:
  file("core_genome.tree") optional true into cgtree

  script:
    """
    numGenomes=`grep -o '>' core_gene_alignment.aln | wc -l`
    if [ \$numGenomes -gt 3 ]
    then
      iqtree -nt AUTO -s core_gene_alignment.aln -m ${params.cg_tree_model} -bb 1000
      mv core_gene_alignment.aln.contree core_genome.tree
    fi
    """
}

//Kraken Step 1: Run Kraken
process kraken {
  tag "$name"
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*_kraken2_report.txt*"

  input:
  set val(name), file(reads) from read_files_kraken

  output:
  tuple name, file("${name}_kraken2_report.txt") into kraken_files, kraken_multiqc

  script:
  """
  kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${name}_kraken2_report.txt --paired ${reads[0]} ${reads[1]}
  """
}

//Kraken Step 2: Summarize kraken results
process kraken_summary {
  tag "$name"
  publishDir "${params.outdir}",mode:'copy'

  input:
  file(files) from kraken_files.collect()

  output:
  file("kraken_results.txt") into kraken_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  files = glob.glob("*kraken2_report*")

  results = []
  for file in files:
      sample_id = os.path.basename(file).split(".")[0].replace("_kraken2_report","")
      df = []
      dfs = []
      with open(file,"r") as inFile:
          for line in inFile:
              line = line.strip()
              sline = line.split("\t")
              if sline[5] == "unclassified":
                  df.append(sline)
              if sline[3] == "S":
                  df.append(sline)

      pandf = DataFrame(df, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
      pandf['Name'] = pandf['Name'].str.lstrip()
      pandf = pandf.sort_values(by=['Percentage'], ascending=False)
      unclass = pandf[pandf["Name"]=="unclassified"]
      pandf = pandf[pandf["Name"]!="unclassified"]
      pandf = pandf.head(2)
      dfs = pd.concat([unclass,pandf])
      dfs = dfs.assign(Sample=sample_id)
      dfs = dfs[['Sample','Percentage','Name']]
      dfs = dfs.reset_index(drop=True)
      print(dfs)

      unclassified = dfs.iloc[0]['Percentage'] + "%"
      first_species = dfs.iloc[1]['Name'] + " (" + dfs.iloc[1]['Percentage'] + "%)"
      second_species = dfs.iloc[2]['Name'] + " (" + dfs.iloc[2]['Percentage'] + "%)"
      combined = [[sample_id, unclassified, first_species, second_species]]
      result = DataFrame(combined, columns=['Sample','Unclassified Reads (%)','Primary Species (%)','Secondary Species (%)'])
      results.append(result)

  df_concat = pd.concat(results)
  df_concat.to_csv(f'kraken_results.txt',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

Channel
  .from("$baseDir/multiqc_config.yaml")
  .set { multiqc_config }

Channel
  .from("$baseDir/assets/dryad_logo_250.png")
  .set { logo }

process multiqc {
  publishDir "${params.outdir}",mode:'copy'

  input:
  file(a) from multiqc_clean_reads.collect()
  file(c) from stats_multiqc.collect()
  file(e) from kraken_multiqc.collect()
  file(f) from logo
  file(config) from multiqc_config
  file(d) from prokka_multiqc.collect()
  file(r) from fastqc_results.collect()

  output:
  file("*.html") into multiqc_output

  script:
  """
  multiqc -c ${config} .
  """
}
