#!/usr/bin/env nextflow

//Description: Workflow for generating a genomic comparison between HAI/AR samples.
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

if(params.test){

  testIDS = ['SRR14311557','SRR14311556','SRR14311555','SRR14311554',
    'SRR14311553','SRR14311552','SRR14613509']

  println "Running test analysis using the following samples:"
  println testIDS

  Channel
      .fromSRA(testIDS)
      .into { raw_reads; raw_reads_count }
  Channel
      .fromPath("$baseDir/assets/ASM211692v1.fasta")
      .into { snp_reference; mapping_reference }
  Channel
      .fromPath("$baseDir/configs/snppipeline.conf")
      .set { snp_config }

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

  if (params.snp_reference) {
      Channel
          .fromPath(params.snp_reference)
          .into { snp_reference; mapping_reference }
  }

  if (params.cfsan_config != "" & params.test != true) {
      Channel
          .fromPath(params.cfsan_config)
          .set { snp_config }
  }

  if (params.cfsan_config == "" & params.test != true) {
      Channel
          .fromPath("$baseDir/configs/snppipeline.conf")
          .set { snp_config }
  }
}

if (params.kraken_db != "") {
    Channel
        .fromPath(params.kraken_db)
        .set { kraken_db }
} else {
    kraken_db = file('NO_FILE')
}

//Preprocessing Step: Change read names
process preProcess {
  //errorStrategy 'ignore'

  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(outfiles) into processed_reads_fastqc, processed_reads_trimming

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

//QC Step: Trim reads and remove PhiX contamination
process clean_reads {
  tag "$name"
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/trimming/stats", mode: 'copy', pattern:"*.trim.txt"
  publishDir "${params.outdir}/trimming/reads", mode: 'copy', pattern:"*.gz"

  input:
  set val(name), file(reads) from processed_reads_trimming

  output:
  tuple name, file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_shovill, cleaned_reads_fastqc, cleaned_reads_mapping, cleaned_reads_kraken
  file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_snp
  file("${name}.adapter.stats.txt") into multiqc_clean_reads
  file("${name}.trim.txt") into bbduk_files

  script:
  """
  bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${name}_clean_1.fastq.gz out2=${name}_clean_2.fastq.gz outm=${name}.adapters.fq qtrim=${params.trimdirection} trimq=${params.qualitytrimscore} minlength=${params.minlength} ref=/bbmap/resources/adapters.fa stats=${name}.adapter.stats.txt k=31 hdist=1 tpe tbo &> ${name}.out
  grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${name}.out > ${name}.trim.txt
  """
}

//QC Step: Summarize BBduk results
process bbduk_summary {
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/trimming", mode:'copy'

  input:
  file(files) from bbduk_files.collect()

  output:
  file("bbduk_results.tsv") into bbduk_tsv

  script:
  """
  #!/usr/bin/python3.7

  import os
  import glob
  import numpy
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
  files = glob.glob("*.trim.txt")

  # summarize bbduk output files
  results = map(summarize_bbduk,files)

  # convert results to data frame and write to tsv
  df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])
  df.to_csv(f'bbduk_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

combined_reads = processed_reads_fastqc.concat(cleaned_reads_fastqc)

//QC Step: Run FastQC
process fastqc {
  tag "$name"
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from combined_reads

  output:
  file("*_fastqc.{zip,html}") into fastqc_results, fastqc_multiqc

  script:
  """
  fastqc -q  ${reads}
  """
}

//QC Step: Summarize FastQC results
process fastqc_summary {
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  file(fastqc) from fastqc_results.collect()

  output:
  file("fastqc_summary.txt") into fastqc_summary

  shell:
  """
  zips=`ls *.zip`
  for i in \$zips; do
      unzip -o \$i &>/dev/null;
  done
  fq_folders=\${zips}
  for folder in \$fq_folders; do
    folder=\${folder%.*}
    cat \$folder/summary.txt >> fastqc_summary.txt
    ls .
  done;
  sed -i 's/.fastq.gz//g' fastqc_summary.txt
  """
}

//Kraken Step 1: Run Kraken
process kraken {
  tag "$name"
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*.kraken2.txt*"

  input:
  set val(name), file(cleaned_reads) from cleaned_reads_kraken
  file(db) from kraken_db

  output:
  tuple name, file("${name}.kraken2.txt") into kraken_files, kraken_multiqc
  file("Kraken2_DB.txt") into kraken_version

  script:
  if (params.kraken_db != "") {
      """
      dbname=${db}
      dbname=\${dbname%.*.*}

      mkdir kraken2-db
      tar -xvf ${db} --directory kraken2-db
      echo \$dbname > Kraken2_DB.txt

      kraken2 --db ./kraken2-db --threads ${task.cpus} --report ${name}.kraken2.txt --paired ${cleaned_reads[0]} ${cleaned_reads[1]}
      """
  } else {
      """
      kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${name}.kraken2.txt --paired ${cleaned_reads[0]} ${cleaned_reads[1]}
      ls /kraken2-db/ > Kraken2_DB.txt
      """
  }
}

//Kraken Step 2: Summarize kraken results
process kraken_summary {
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/kraken", mode:'copy'

  input:
  file(files) from kraken_files.collect()

  output:
  file("kraken_results.tsv") into kraken_tsv
  file("kraken_results.tsv") into kraken_prokka

  script:
  """
  #!/usr/bin/python3.7

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
  files = glob.glob("*.kraken2.txt*")

  # summarize kraken2 report files
  results = map(summarize_kraken, files)

  # concatenate summary results and write to tsv
  data_concat = pd.concat(results)
  data_concat.to_csv(f'kraken_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//CG Step: Assemble cleaned reads with Shovill and map reads back to the assembly with BWA
process shovill {
  tag "$name"
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/assembled", mode: 'copy', pattern:"*.fa"
  publishDir "${params.outdir}/mapping/sams", mode: 'copy', pattern:"*.sam"

  input:
  set val(name), file(cleaned_reads) from cleaned_reads_shovill

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_prokka
  tuple name, file("${name}.assembly.sam") into assembly_sams

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
  tag "$name"
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/quast", mode:'copy', pattern: "*.quast.report.tsv"

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.transposed.quast.report.tsv") into quast_files
  file("${name}.quast.report.tsv") into quast_multiqc

  script:
  """
  quast.py ${name}.contigs.fa -o .
  mv report.tsv ${name}.quast.report.tsv
  mv transposed_report.tsv ${name}.transposed.quast.report.tsv
  """
}

//QC Step: Summarize QUAST results
process quast_summary {
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/quast", mode:'copy'

  input:
  file(files) from quast_files.collect()

  output:
  file("quast_results.tsv") into quast_tsv

  script:
  """
  #!/usr/bin/python3.7

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
  files = glob.glob("*.transposed.quast.report.tsv")

  # summarize quast output files
  dfs = map(summarize_quast,files)
  dfs = list(dfs)

  # concatenate dfs and write data frame to file
  dfs_concat = pd.concat(dfs)
  dfs_concat.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//CG Step: Set up Prokka files
process prokka_setup {
  tag "$name"
  //errorStrategy 'ignore'

  input:
  file(kraken) from kraken_prokka
  set val(name), file(input) from assembled_genomes_prokka

  output:
  tuple name, file("${name}.*.fa") into prokka_input

  script:
  """
  #!/usr/bin/python3.7

  import os
  import pandas as pd
  import shutil

  genomeFile = '${input}'
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

  publishDir "${params.outdir}/annotated", mode:'copy'

  input:
  set val(name), file(assembly) from prokka_input

  output:
  file("${name}.gff") into annotated_genomes
  file("${name}.prokka.stats.txt") into prokka_multiqc

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
  //errorStrategy 'ignore'

  publishDir "${params.outdir}", mode:'copy'

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

//CG Step: Infer ML tree from core genome alignment using IQ-TREE
process cg_tree {
  //errorStrategy 'ignore'

  publishDir "${params.outdir}", mode:'copy'

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

//Run SNP calling pipeline if a reference sequence is provided
if (params.snp_reference != null & !params.snp_reference.isEmpty() | params.test_snp & params.test) {

    //SNP Step: Run CFSAN-SNP Pipeline
    process cfsan {
      //errorStrategy 'ignore'

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
      //errorStrategy 'ignore'

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

    //SNP Step: Map cleaned reads to reference sequence using BWA
    process bwa {
      tag "$name"
      //errorStrategy 'ignore'

      publishDir "${params.outdir}/mapping/sams", mode: 'copy', pattern:"*.sam"

      input:
      file(reference) from mapping_reference.first()
      set val(name), file(reads) from cleaned_reads_mapping

      output:
      tuple name, file("${name}.reference.sam") into reference_sams

      script:
      """
      bwa index ${reference}
      bwa mem ${reference} ${reads[0]} ${reads[1]} > ${name}.reference.sam
      """
    }
    //Combine assembly and reference alignment channels
    sam_files = assembly_sams.concat(reference_sams)
}

else {
  //Otherwise use only assembly alignment files
  sam_files = assembly_sams
}

//QC Step: Convert SAMs to BAMs, index and sort BAM files, then calculate coverage
process samtools {
  tag "$name"
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/mapping/bams", mode: 'copy', pattern:"*.sorted.bam*"
  publishDir "${params.outdir}/mapping/depth", mode: 'copy', pattern:"*.depth.tsv"
  publishDir "${params.outdir}/mapping/stats", mode: 'copy', pattern:"*.stats.txt"

  input:
  set val(name), file(sam) from sam_files

  output:
  tuple name, file("*.depth.tsv") into reference_depth_results,assembly_depth_results
  tuple name, file("*.mapped.tsv") into reference_mapped_results,assembly_mapped_results
  file("*.stats.txt")
  file("*.sorted.bam*")

  shell:
  """
  filename=${sam}
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
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/mapping", mode: 'copy'

  input:
  file(depth) from assembly_depth_results.collect()
  file(mapped) from assembly_mapped_results.collect()

  output:
  file('coverage_stats.tsv') into assembly_mapping_tsv

  script:
  """
  #!/usr/bin/python3.7

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
  files = glob.glob("*.assembly.depth.tsv*")

  # summarize samtools depth files
  results = map(summarize_depth,files)

  # write results to file
  with open('coverage_stats.tsv', 'w') as outFile:
      outFile.write("Sample\\tMedian Coverage (Mapped to Assembly)\\tMean Coverage (Mapped to Assembly)\\n")
      for result in results:
          outFile.write(result)
  """
}

if (params.snp_reference != null & !params.snp_reference.isEmpty() | params.test_snp & params.test) {

    //QC Step: Calculate coverage for reads mapped to reference
    process reference_mapping_stats {
      //errorStrategy 'ignore'

      publishDir "${params.outdir}/mapping", mode: 'copy'

      input:
      file(depth) from reference_depth_results.collect()
      file(mapped) from reference_mapped_results.collect()

      output:
      file('mapping_results.tsv') into reference_mapping_tsv

      script:
      """
      #!/usr/bin/python3.7

      import pandas as pd
      import os
      import glob
      from functools import reduce

      depth_files = glob.glob("*.reference.depth.tsv")
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

      read_files = glob.glob("*.reference.mapped.tsv")
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
}
else {
    reference_mapping_tsv = Channel.empty()
}

//QC Step: Merge QC results into one tsv
process merge_results {
  //errorStrategy 'ignore'

  publishDir "${params.outdir}/", mode: 'copy'

  input:
  file(bbduk) from bbduk_tsv
  file(quast) from quast_tsv
  file(assembly) from assembly_mapping_tsv
  file(kraken) from kraken_tsv
  file(vkraken) from kraken_version.first()
  file(reference) from reference_mapping_tsv.ifEmpty{ 'empty' }

  output:
  file('dryad_report.csv')

  script:
  """
  #!/usr/bin/python3.7

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
  .fromPath("$baseDir/configs/multiqc_config.yaml")
  .set { multiqc_config }

Channel
  .fromPath("$baseDir/assets/dryad_logo_250.png")
  .set { logo }

process multiqc {
  //errorStrategy 'ignore'

  publishDir "${params.outdir}", mode:'copy'

  input:
  file(a) from multiqc_clean_reads.collect()
  file(b) from fastqc_multiqc.collect()
  file(c) from kraken_multiqc.collect()
  file(d) from quast_multiqc.collect()
  file(e) from prokka_multiqc.collect()
  file(f) from logo
  file(config) from multiqc_config

  output:
  file("*.html") into multiqc_output

  script:
  """
  multiqc -c ${config} .
  """
}
