#!/usr/bin/env nextflow

//Description: Workflow for generating a genomic comparison between HAI/AR samples.
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

nextflow.enable.dsl=2

params.test = false
params.test_snp = true
if(params.test){
  testIDS = ['SRR14311557','SRR14311556','SRR14311555','SRR14311554',
    'SRR14311553','SRR14311552','SRR14613509']
  println "Running test analysis using the following samples:"
  println testIDS
  Channel
      .fromSRA(testIDS)
      .set { raw_reads }
  Channel
      .fromPath("$baseDir/assets/ASM211692v1.fasta")
      .set { snp_reference }
  Channel
      .fromPath("$baseDir/configs/snppipeline.conf")
      .set { snp_config }
} else{
  //setup channel to read in and pair the fastq files
  Channel
      .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
      .set { raw_reads }
  //check we have at least 3 samples
  Channel
      .from(raw_reads)
      .collect()
      .subscribe {
        int size = it.queue[0].size()
        if(size < 3){
          println "Dryad requires 3 or more samples."
          System.exit(1)
        }
      }
  if (params.snp_reference) {
      Channel
          .fromPath(params.snp_reference)
          .set { snp_reference }
      Channel
          .fromPath("$baseDir/configs/snppipeline.conf")
          .set { snp_config }
  }
}

//Preprocessing Step: Change read names
process preProcess {
  //publishDir "${params.outdir}/reads", mode: 'copy', pattern:"*.gz"

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("*_R{1,2}.fastq.gz"), emit: processed_reads

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

//QC Step: Trim reads and remove adapters and remove PhiX contamination
process clean_reads {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/trimming/stats", mode: 'copy', pattern:"*.trim.txt"
  publishDir "${params.outdir}/trimming/reads", mode: 'copy', pattern:"*.gz"

  input:
  tuple val(name), path(processed_reads)

  output:
  tuple val(name), path("${name}_clean{_1,_2}.fastq.gz"), emit: cleaned_reads
  path("${name}.trim.txt"), emit: bbduk_files
  path("${name}.adapter.stats.txt"), emit: bbduk_stats
  path("${name}_clean{_1,_2}.fastq.gz"), emit: cleaned_reads_cfsan

  script:
  """
  bbduk.sh in1=${processed_reads[0]} in2=${processed_reads[1]} out1=${name}_clean_1.fastq.gz out2=${name}_clean_2.fastq.gz outm=${name}.adapters.fq qtrim=${params.trimdirection} trimq=${params.qualitytrimscore} minlength=${params.minlength} ref=/bbmap/resources/adapters.fa stats=${name}.adapter.stats.txt k=31 hdist=1 tpe tbo &> ${name}.out
  grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${name}.out > ${name}.trim.txt
  """
}

//Summary Step: Summarize BBDuk results
process bbduk_summary {
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/trimming",mode:'copy'

  input:
  path("data*/*")

  output:
  path("bbduk_results.tsv"), emit: bbduk_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  # function for summarizing bbduk output
  def summarize_bbduk(file):
      # get sample id from file name and set up data list
      sample_id = os.path.basename(file).split(".")[0]
      data = []
      data.append(sample_id)
      with open(file,"r") as inFile:
          for i, line in enumerate(inFile):
              # get total number of reads
              if i == 0:
                  num_reads = line.strip().split("\\t")[1].replace(" reads ","")
                  data.append(num_reads)
              # get total number of reads removed
              if i == 3:
                  rm_reads = line.strip().split("\\t")[1].replace("reads ","")
                  rm_reads = rm_reads.rstrip()
                  data.append(rm_reads)
      return data

  # get all bbduk output files
  files = glob.glob("data*/*.trim.txt")

  # summarize bbduk output files
  results = map(summarize_bbduk,files)

  # convert results to data frame and write to tsv
  df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])
  df.to_csv(f'bbduk_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//QC Step: Run FastQC on processed and cleaned reads
process fastqc {
  //errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  tuple val(name), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//Summary Step: Summarize FastQC results
process fastqc_summary {
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  path(fastqc_results)

  output:
  path("fastqc_summary.tsv"), emit: fastqc_summary

  shell:
  """
  zips=`ls *.zip`

  for i in \$zips; do
      unzip -o \$i &>/dev/null;
  done

  fq_folders=\${zips}

  for folder in \$fq_folders; do
    folder=\${folder%.*}
    cat \$folder/summary.txt >> fastqc_summary.tsv
    ls .
  done;

  sed -i 's/.fastq.gz//g' fastqc_summary.tsv
  """
}

//Classification Step: Run Kraken
process kraken {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*.kraken2.txt*"

  input:
  tuple val(name), path(cleaned_reads)

  output:
  path("${name}.kraken2.txt"), emit: kraken_reports
  path("Kraken2_DB.txt"), emit: kraken_version

  script:
  """
  kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${name}.kraken2.txt --paired ${cleaned_reads[0]} ${cleaned_reads[1]}

  ls /kraken2-db/ > Kraken2_DB.txt
  """
}

//Summary Step: Summarize kraken results
process kraken_summary {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/kraken",mode:'copy'

  input:
  path("data*/*")

  output:
  path("kraken_results.tsv"), emit: kraken_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  # function for summarizing kraken2 report files
  def summarize_kraken(file):
      # get sample id from file name
      sample_id = os.path.basename(file).split('.')[0].replace('.kraken2.txt','')
      data = []
      # read kraken2 report file
      with open(file,'r') as inFile:
          for line in inFile:
              line = line.strip()
              sline = line.split('\\t')
              # get unclassified reads result (denoted by 'unclassified') and append to data
              if sline[5] == 'unclassified':
                  data.append(sline)
              # get species results (denoted by 'S') and append to data
              if sline[3] == 'S':
                  data.append(sline)
      # convert data list to data frame
      data_df = DataFrame(data, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
      # remove left leading spaces from the Name column
      data_df['Name'] = data_df['Name'].str.lstrip()
      # sort data frame by percentages (largest to smallest)
      data_df = data_df.sort_values(by=['Percentage'], ascending=False)
      # make new data frame for unclassified reads only
      unclass = data_df[data_df['Name']=='unclassified']
      # exception for if no unclassified reads found
      if unclass.empty:
          # import pandas as pd
          lst = [['0','NA','NA','NA','NA','NA']]
          unclass = pd.DataFrame(lst, columns =['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
      # subset data frame by species
      species_df = data_df[data_df['Name']!='unclassified']
      # get first two species matches (first two largest percentages) in data frame
      species_df = species_df.head(2)
      # check if species data frame has two rows
      if len(species_df) == 0:
          # add two empty rows to species data frame
          species_df = species_df.append(pd.Series(), ignore_index=True)
          species_df = species_df.append(pd.Series(), ignore_index=True)
      if len(species_df) == 1:
          # add one empty row to species data frame
          species_df = species_df.append(pd.Series(), ignore_index=True)
      # concatenate unclassified data frame and species data frame
      df_concat = pd.concat([unclass,species_df])
      # add sample name column to concatenated data frame
      df_concat = df_concat.assign(Sample=sample_id)
      # keep only Sample Percentage and Name columns in concatenated data frame
      df_concat = df_concat[['Sample','Percentage','Name']]
      # reset index of concatenated data frame using drop parameter to avoid old index added as column
      df_concat = df_concat.reset_index(drop=True)
      # add percentage sign to unclassified column
      unclassified = df_concat.iloc[0]['Percentage'] + '%'
      # convert to lists
      # if primary species is nan, replace with NA
      if str(df_concat.iloc[1]['Name']) == 'nan':
          primary_species = 'NA'
      # otherwise convert to (#%)
      else:
          primary_species = df_concat.iloc[1]['Name'] + ' (' + df_concat.iloc[1]['Percentage'] + '%)'
      # repeat for secondary species
      if str(df_concat.iloc[2]['Name']) == 'nan':
          secondary_species = 'NA'
      else:
          secondary_species = df_concat.iloc[2]['Name'] + ' (' + df_concat.iloc[2]['Percentage'] + '%)'
      # list of lists
      combined = [[sample_id, unclassified, primary_species, secondary_species]]
      # convert list of lists to data frame
      combined_df = DataFrame(combined, columns=['Sample','Unclassified Reads (%)','Primary Species (%)','Secondary Species (%)'])
      return combined_df
  # get all kraken2 report files
  files = glob.glob("data*/*.kraken2.txt")
  # summarize kraken2 report files
  results = map(summarize_kraken, files)
  # concatenate summary results and write to tsv
  data_concat = pd.concat(results)
  data_concat.to_csv(f'kraken_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}


//Assembly step: Assemble trimmed reads with Shovill and map reads back to assembly
process shovill {
  //errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/assembled", mode: 'copy', pattern:"*.fa"
  publishDir "${params.outdir}/mapping/sams", mode: 'copy', pattern:"*.sam"

  input:
  tuple val(name), path(cleaned_reads)

  output:
  tuple val(name), path("${name}.contigs.fa"), emit: assembled_genomes
  tuple val(name), path("${name}.assembly.sam"), emit: assembly_sams

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory} --outdir ./output --R1 ${cleaned_reads[0]} --R2 ${cleaned_reads[1]} --force
  mv ./output/contigs.fa ${name}.contigs.fa
  bwa index ${name}.contigs.fa
  bwa mem ${name}.contigs.fa ${cleaned_reads[0]} ${cleaned_reads[1]} > ${name}.assembly.sam
  """
}

//QC Step: Run QUAST on assemblies
process quast {
  //errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/quast",mode:'copy', pattern: "${name}.transposed.quast.tsv"

  input:
  tuple val(name), path(assembled_genomes)

  output:
  path("*.transposed.quast.tsv"), emit: quast_files
  path("*.report.quast.tsv"), emit: quast_reports

  script:
  """
  quast.py ${name}.contigs.fa -o .
  mv report.tsv ${name}.report.quast.tsv
  mv transposed_report.tsv ${name}.transposed.quast.tsv
  """
}

//Summary Step: Summarize QUAST results
process quast_summary {
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  path("data*/*")

  output:
  path("quast_results.tsv"), emit: quast_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  # function for summarizing quast output
  def summarize_quast(file):
      # get sample id from file name and set up data list
      sample_id = os.path.basename(file).split(".")[0]
      # read in data frame from file
      df = pd.read_csv(file, sep='\\t')
      # get contigs, total length and assembly length columns
      df = df.iloc[:,[1,7,17]]
      # assign sample id as column
      df = df.assign(Sample=sample_id)
      # rename columns
      df = df.rename(columns={'# contigs (>= 0 bp)':'Contigs','Total length (>= 0 bp)':'Assembly Length (bp)'})
      # re-order data frame
      df = df[['Sample', 'Contigs','Assembly Length (bp)', 'N50']]
      return df

  # get quast output files
  files = glob.glob("data*/*.transposed.quast.tsv")

  # summarize quast output files
  dfs = map(summarize_quast,files)

  # concatenate dfs and write data frame to file
  dfs_concat = pd.concat(dfs)
  dfs_concat.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//CG Step: Set up files for Prokka
process prokka_setup {
  tag "$name"

  input:
  path("kraken_results.tsv")
  tuple val(name), path(assembled_genome)

  output:
  tuple val(name), path("${name}.*.fa"), emit: prokka_input

  script:
  """
  #!/usr/bin/env python3
  import os
  import pandas as pd
  import shutil

  genomeFile = '${assembled_genome}'
  sid = genomeFile.split('.')[0]
  df = pd.read_csv('kraken_results.tsv', header=0, delimiter='\\t')
  df = df[df['Sample'] == sid]
  taxa = df.iloc[0]['Primary Species (%)']
  taxa = taxa.split(' ')
  taxa = taxa[0] + '_' + taxa[1]
  shutil.copyfile(genomeFile, f'{sid}.{taxa}.fa')
  """
}

//CG Step: Run Prokka on assemblies
process prokka {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/annotated",mode:'copy'

  input:
  tuple val(name), path(assembly)

  output:
  path("${name}.gff"), emit: annotated_genomes
  path("${name}.prokka.stats.txt"), emit: prokka_stats

  script:
  """
  filename=${assembly}
  handle=\${filename%.*}
  taxa=\${handle##*.}
  genus=\${taxa%_*}
  species=\${taxa##*_}

  prokka --cpu ${task.cpus} --force --compliant --prefix ${name} --genus \$genus --species \$species --strain ${name} --mincontiglen 500 --outdir . ${assembly} > ${name}.log
  mv ${name}.txt ${name}.prokka.stats.txt
  """
}

//CG Step: Perform core genome alignment using Roary
process roary {
  publishDir "${params.outdir}",mode:'copy'

  numGenomes = 0
  input:
  path("data*/*")

  output:
  path("core_gene_alignment.aln"), emit: core_genome_alignment
  path("core_genome_statistics.txt"), emit: core_genome_statistics

  script:
  if(params.roary_mafft == true){
    mafft="-n"
  }else{mafft=""}
  """
  cpus=`grep -c ^processor /proc/cpuinfo`
  roary -e ${mafft} -p \$cpus data*/*.gff
  mv summary_statistics.txt core_genome_statistics.txt
  """
}

//CG Step: Infer ML tree from core genome alignment using IQ-TREE
process cg_tree {
  publishDir "${params.outdir}",mode:'copy'

  input:
  path(alignment)

  output:
  path("core_genome.tree"), optional: true, emit: cgtree

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

//Run SNP calling pipeline if a reference sequence is provided
// if (params.snp_reference != null & !params.snp_reference.isEmpty() | params.test_snp & params.test) {

//SNP Step: Run CFSAN-SNP Pipeline
process cfsan {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path("data*/*")
  path(reference)
  path(config)

  output:
  path("snp_distance_matrix.tsv"), emit: snp_matrix
  path("snpma.fasta"), emit: snp_alignment

  script:
  """
  #!/usr/bin/env python
  import subprocess
  import glob
  import os,sys
  import shutil

  fwd_reads = glob.glob("data*/*_clean_1.fastq.gz*")
  fwd_reads.sort()
  rev_reads = glob.glob("data*/*_clean_2.fastq.gz*")
  rev_reads.sort()

  if len(fwd_reads) != len(rev_reads):
    sys.exit("Uneven number of forward and reverse reads.")

  os.mkdir("input_reads")

  c = 0
  while c < len(fwd_reads):
    name = os.path.basename(fwd_reads[c]).split('_clean_1')[0]
    path = os.path.join("input_reads",name)
    os.mkdir(path)
    new_fwd_path = os.path.join(path,os.path.basename(fwd_reads[c]))
    new_rev_path = os.path.join(path,os.path.basename(rev_reads[c]))
    shutil.copy(fwd_reads[c],new_fwd_path)
    shutil.copy(rev_reads[c],new_rev_path)
    #os.rename(fwd_reads[c],new_fwd_path)
    #os.rename(rev_reads[c],new_rev_path)
    c += 1

    # name = os.path.basename(fwd_reads[c]).split('_clean_1')[0]
    # path = os.path.join("input_reads",name)
    # os.mkdir(path)
    # os.rename(fwd_reads[c],os.path.join(path,fwd_reads[c]))
    # os.rename(rev_reads[c],os.path.join(path,rev_reads[c]))
    # c += 1

  command = "cfsan_snp_pipeline run ${reference} -c ${config} -o . -s input_reads"
  # command = "cfsan_snp_pipeline run ${reference} -o . -s input_reads"
  process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
  output, error = process.communicate()
  print output
  print error
  """
}

//SNP Step: Run IQ-TREE on SNP alignment
process snp_tree {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(snp_alignment)

  output:
  path("snp.tree"), optional: true, emit: snptree

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

//SNP Step: Map cleaned reads to reference sequence using BWA
process bwa {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/mapping/sams", mode: 'copy', pattern:"*.sam"

  input:
  path(reference)
  tuple val(name), path(cleaned_reads)

  output:
  tuple val(name), path("${name}.reference.sam"), emit: reference_sams

  script:
  """
  bwa index ${reference}
  bwa mem ${reference} ${cleaned_reads[0]} ${cleaned_reads[1]} > ${name}.reference.sam
  ls .
  """
}

//QC Step: Convert SAMs to BAMs, index and sort BAM files, then calculate coverage
process samtools {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/mapping/bams", mode: 'copy', pattern:"*.sorted.bam*"
  publishDir "${params.outdir}/mapping/depth", mode: 'copy', pattern:"*.depth.tsv"
  publishDir "${params.outdir}/mapping/stats", mode: 'copy', pattern:"*.stats.txt"

  input:
  tuple val(name), path(sam_files)

  output:
  path("*.depth.tsv"), emit: depth_results
  path("*.mapped.tsv"), emit: mapping_results
  path("*.stats.txt"), emit: samtools_stats
  path("*.sorted.bam*")

  shell:
  """
  filename=${sam_files}
  handle=\${filename%.*}
  type=\${handle##*.}

  samtools view -S -b ${name}.\$type.sam > ${name}.\$type.bam
  samtools sort ${name}.\$type.bam > ${name}.\$type.sorted.bam
  samtools index ${name}.\$type.sorted.bam
  samtools depth -a ${name}.\$type.sorted.bam > ${name}.\$type.depth.tsv
  samtools stats ${name}.\$type.sorted.bam > ${name}.\$type.stats.txt
  samtools view -c -F 260 ${name}.\$type.sorted.bam > ${name}.\$type.mapped.tsv
  samtools view -c ${name}.\$type.sorted.bam >> ${name}.\$type.mapped.tsv
  """
}

//QC Step: Using samtools depth, calculate coverage of cleaned reads mapped to their assemblies
process assembly_coverage_stats {
  publishDir "${params.outdir}/mapping", mode: 'copy'

  input:
  path("data*/*")

  output:
  path('coverage_stats.tsv'), emit: assembly_mapping_tsv

  script:
  """
  #!/usr/bin/env python3
  import glob
  import os
  from numpy import median
  from numpy import average

  # function for summarizing samtools depth files
  def summarize_depth(file):
      # get sample id from file name and set up data list
      sid = os.path.basename(file).split('.')[0]
      data = []
      # open samtools depth file and get depth
      with open(file,'r') as inFile:
          for line in inFile:
              data.append(int(line.strip().split()[2]))
      # get median and average depth
      med = int(median(data))
      avg = int(average(data))
      # return sample id, median and average depth
      result = f"{sid}\\t{med}\\t{avg}\\n"
      return result

  # get all samtools depth files
  files = glob.glob("data*/*.assembly.depth.tsv*")

  # summarize samtools depth files
  results = map(summarize_depth,files)

  # write results to file
  with open('coverage_stats.tsv', 'w') as outFile:
      outFile.write("Sample\\tMedian Coverage (Mapped to Assembly)\\tMean Coverage (Mapped to Assembly)\\n")
      for result in results:
          outFile.write(result)
  """
}

//QC Step: Calculate coverage for reads mapped to reference
process reference_mapping_stats {
  publishDir "${params.outdir}/mapping", mode: 'copy'

  input:
  path("data*/*")

  output:
  path("mapping_results.tsv"), emit: reference_mapping_tsv

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import os
  import glob
  from functools import reduce

  depth_files = glob.glob("data*/*.reference.depth.tsv")
  depth_dfs = []
  cols = ["Sample","Base Pairs Mapped to Reference >1X (%)","Base Pairs Mapped to Reference >40X (%)"]
  depth_dfs.append(cols)

  for file in depth_files:
      sampleID = os.path.basename(file).split(".")[0]
      depth_df = pd.read_csv(file, sep="\\t", header=None)
      overForty = int((len(depth_df[(depth_df[2]>40)])/len(depth_df)) * 100)
      overOne = int((len(depth_df[(depth_df[2]>1)])/len(depth_df)) * 100)
      stats = [sampleID, overOne, overForty]
      depth_dfs.append(stats)
  depth_df = pd.DataFrame(depth_dfs[1:], columns=depth_dfs[0])

  read_files = glob.glob("data*/*.reference.mapped.tsv")
  read_dfs = []
  cols = ["Sample","Reads Mapped to Reference (%)"]
  read_dfs.append(cols)
  for file in read_files:
      sampleID = os.path.basename(file).split(".")[0]
      read_df = pd.read_csv(file, sep="\\t", header=None)
      mapped_reads = read_df.iloc[0][0]
      all_reads = read_df.iloc[1][0]
      percent_mapped = int((mapped_reads/all_reads) * 100)
      stats = [sampleID,percent_mapped]
      read_dfs.append(stats)
  read_dfs = pd.DataFrame(read_dfs[1:], columns=read_dfs[0])

  dfs = [depth_df, read_dfs]
  merged = reduce(lambda  left,right: pd.merge(left,right,on=["Sample"], how="left"), dfs)
  merged[['Reads Mapped to Reference (%)','Base Pairs Mapped to Reference >1X (%)','Base Pairs Mapped to Reference >40X (%)']] = merged[['Reads Mapped to Reference (%)','Base Pairs Mapped to Reference >1X (%)','Base Pairs Mapped to Reference >40X (%)']].astype(str) + '%'
  merged.to_csv("mapping_results.tsv",sep="\\t", index=False, header=True, na_rep="NaN")
  """
}

//QC Step: Merge QC results into one tsv
process merge_results {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  path("bbduk_results.tsv")
  path("quast_results.tsv")
  path("coverage_stats.tsv")
  path("kraken_results.tsv")
  path("Kraken2_DB.txt")
  path("mapping_results.tsv")

  output:
  file('dryad_report.csv')

  script:
  """
  #!/usr/bin/env python3

  import os
  import glob
  import pandas as pd
  from functools import reduce

  with open('Kraken2_DB.txt', 'r') as krakenFile:
      krakenDB_version = krakenFile.readline().strip()


  files = glob.glob('*.tsv')

  dfs = []

  for file in files:
      df = pd.read_csv(file, header=0, delimiter='\\t')
      dfs.append(df)

  merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],
                                              how='left'), dfs)
  merged = merged.assign(krakenDB=krakenDB_version)

  merged = merged.rename(columns={'Contigs':'Contigs (#)','krakenDB':'Kraken Database Verion'})

  merged.to_csv('dryad_report.csv', index=False, sep=',', encoding='utf-8')
  """
}

Channel
  .fromPath("$baseDir/multiqc_config.yaml")
  .set { multiqc_config }

Channel
  .fromPath("$baseDir/assets/dryad_logo_250.png")
  .set { logo }

//Summary Step: MultiQC
process multiqc {
  publishDir "${params.outdir}",mode:'copy'

  input:
  path("data*/*")
  path(config)
  path(logo)

  output:
  path("*.html"), emit: multiqc_output

  script:
  """
  multiqc -c ${config} .
  """
}

workflow {

    preProcess(raw_reads)

    clean_reads(preProcess.out.processed_reads)

    bbduk_summary(clean_reads.out.bbduk_files.collect())

    processed = preProcess.out.processed_reads
    cleaned = clean_reads.out.cleaned_reads
    combined_reads = processed.concat(cleaned)

    fastqc(combined_reads)

    fastqc_summary(fastqc.out.fastqc_results.collect())

    shovill(clean_reads.out.cleaned_reads)

    quast(shovill.out.assembled_genomes)

    quast_summary(quast.out.quast_files.collect())

    kraken(clean_reads.out.cleaned_reads)

    kraken_summary(kraken.out.kraken_reports.collect())

    prokka_setup(kraken_summary.out.kraken_tsv,shovill.out.assembled_genomes)

    prokka(prokka_setup.out.prokka_input)

    roary(prokka.out.annotated_genomes.collect())

    cg_tree(roary.out.core_genome_alignment)

    if (params.snp_reference != null & !params.snp_reference.isEmpty() | params.test_snp & params.test) {
      cfsan(clean_reads.out.cleaned_reads_cfsan.collect(),snp_reference,snp_config)

      snp_tree(cfsan.out.snp_alignment)

      bwa(snp_reference.first(),clean_reads.out.cleaned_reads)

      assembly_sams = shovill.out.assembly_sams
      reference_sams = bwa.out.reference_sams
      sam_files = assembly_sams.concat(reference_sams)
    }
    else {
      sam_files = shovill.out.assembly_sams
    }

    samtools(sam_files)
    assembly_coverage_stats(samtools.out.depth_results.collect())

    if (params.snp_reference != null & !params.snp_reference.isEmpty() | params.test_snp & params.test) {
      depth = samtools.out.depth_results
      mapping = samtools.out.mapping_results
      reference_samtools_results = depth.concat(mapping)
      reference_mapping_stats(reference_samtools_results.collect())
      reference_mapping_tsv = reference_mapping_stats.out.reference_mapping_tsv
    }
    else {
      reference_mapping_tsv = Channel.empty()
    }

    merge_results(bbduk_summary.out.bbduk_tsv,quast_summary.out.quast_tsv,assembly_coverage_stats.out.assembly_mapping_tsv,kraken_summary.out.kraken_tsv,kraken.out.kraken_version.first(),reference_mapping_tsv.ifEmpty{ 'empty' })

    multiqc(clean_reads.out.bbduk_stats.mix(clean_reads.out.bbduk_stats,fastqc.out.fastqc_results,kraken.out.kraken_reports,quast.out.quast_reports,prokka.out.prokka_stats,samtools.out.samtools_stats).collect(),multiqc_config,logo)
}
