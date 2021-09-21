#!/usr/bin/env nextflow

//Description: Workflow for generating a genomic comparison between HAI/AR samples.
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

params.test = false
params.test_snp = true
if(params.test){
  testIDS = ['SRR14311557','SRR14311556','SRR14311555','SRR14311554',
    'SRR14311553','SRR14311552','SRR14613509','SRR14874874']
  println "Running test analysis using the following samples:"
  println testIDS
  Channel
      .fromSRA(testIDS)
      .into { raw_reads; raw_reads_count }
  Channel
      .fromPath("$baseDir/assets/ASM211692v1.fasta")
      .into { snp_reference;mapping_reference }
  Channel
      .from("$baseDir/snppipeline.conf")
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
          .into { snp_reference;mapping_reference }
      Channel
          .from("$baseDir/snppipeline.conf")
          .set { snp_config }
  }
}

//Preprocessing Step: Change read names
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

//QC Step: Trim reads and remove PhiX contamination
process clean_reads {
  tag "$name"
  publishDir "${params.outdir}/trimming", mode: 'copy',pattern:"*.trim.txt"

  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_shovill, cleaned_reads_fastqc, cleaned_reads_mapping, cleaned_reads_kraken
  file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_snp
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats
  tuple file("${name}.phix.stats.txt"),file("${name}.adapters.stats.txt"),file("${name}.trim.txt") into multiqc_clean_reads
  file("${name}.trim.txt") into bbduk_files

  script:
  """
  bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.trimmed_1.fastq.gz out2=${name}.trimmed_2.fastq.gz qtrim=window,${params.windowsize} trimq=${params.qualitytrimscore} minlength=${params.minlength} tbo tbe &> ${name}.out
  repair.sh in1=${name}.trimmed_1.fastq.gz in2=${name}.trimmed_2.fastq.gz out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_clean_1.fastq.gz out2=${name}_clean_2.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${name}.out > ${name}.trim.txt
  """
}

//QC Step: Summarize BBduk results
process bbduk_summary {
  publishDir "${params.outdir}/trimming",mode:'copy'

  input:
  file(files) from bbduk_files.collect()

  output:
  file("bbduk_results.tsv") into bbduk_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  files = glob.glob("*.txt")

  results = []
  for file in files:
      sample_id = os.path.basename(file).split(".")[0]
      vals = []
      vals.append(sample_id)
      with open(file,"r") as inFile:
          for i, line in enumerate(inFile):
              if i == 0:
                  num_reads = line.strip().split("\\t")[1].replace(" reads ","")
                  vals.append(num_reads)
              if i == 3:
                  rm_reads = line.strip().split("\\t")[1].replace("reads ","")
                  rm_reads = rm_reads.rstrip()
                  vals.append(rm_reads)
      results.append(vals)

  df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])
  df.to_csv(f'bbduk_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

combined_reads = read_files_fastqc.concat(cleaned_reads_fastqc)

//QC Step: Run FastQC
process fastqc {
  tag "$name"
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

//Species identification step: Run Kraken
process kraken {
  tag "$name"
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*_kraken2_report.txt*"

  input:
  set val(name), file(reads) from cleaned_reads_kraken

  output:
  tuple name, file("${name}_kraken2_report.txt") into kraken_files, kraken_multiqc

  script:
  """
  kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${name}_kraken2_report.txt --paired ${reads[0]} ${reads[1]}
  """
}

//Species identification step: Summarize Kraken results
process kraken_summary {
  tag "$name"
  publishDir "${params.outdir}/kraken",mode:'copy'

  input:
  file(files) from kraken_files.collect()

  output:
  file("kraken_results.tsv") into kraken_tsv
  file("kraken_results.tsv") into kraken_prokka

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
  df_concat.to_csv(f'kraken_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//CG Step: Assemble cleaned reads with Shovill and map reads back to the assembly with BWA
process shovill {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/assembled", mode: 'copy',pattern:"*.fa"
  publishDir "${params.outdir}/mapping/sams", mode: 'copy',pattern:"*.sam"

  input:
  set val(name), file(reads) from cleaned_reads_shovill

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_prokka
  tuple name, file("${name}.assembly.sam") into assembly_sams

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory} --outdir ./output --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv ./output/contigs.fa ${name}.contigs.fa
  bwa index ${name}.contigs.fa
  bwa mem ${name}.contigs.fa ${reads[0]} ${reads[1]} > ${name}.assembly.sam
  """
}


//QC Step: Run QUAST on assemblies
process quast {
  tag "$name"

  errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy',pattern: "${name}.quast.tsv"

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_files
  file("${name}.report.quast.tsv") into quast_multiqc

  script:
  """
  quast.py ${assembly} -o .
  mv report.tsv ${name}.report.quast.tsv
  mv transposed_report.tsv ${name}.quast.tsv
  """
}

//QC Step: Summarize QUAST results
process quast_summary {
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  file(files) from quast_files.collect()

  output:
  file("quast_results.tsv") into quast_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  files = glob.glob("*.quast.tsv")

  dfs = []

  for file in files:
      sample_id = os.path.basename(file).split(".")[0]
      df = pd.read_csv(file, sep='\\t')
      df = df.iloc[:,[1,7,17]]
      df = df.assign(Sample=sample_id)
      df = df.rename(columns={'# contigs (>= 0 bp)':'Contigs','Total length (>= 0 bp)':'Assembly Length (bp)'})
      df = df[['Sample', 'Contigs','Assembly Length (bp)', 'N50']]
      dfs.append(df)

  dfs_concat = pd.concat(dfs)
  dfs_concat.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//CG Step: Set up Prokka files
process prokka_setup {
  tag "$name"

  input:
  file(kraken) from kraken_prokka
  set val(name), file(input) from assembled_genomes_prokka

  output:
  tuple name, file("${name}.*.fa") into prokka_input

  script:
  """
  #!/usr/bin/env python3
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
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/annotated",mode:'copy'

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

  prokka --cpu ${task.cpus} --force --compliant --prefix ${name} --genus \$genus --species \$species --mincontiglen 500 --outdir . ${assembly} > ${name}.log
  mv ${name}.txt ${name}.prokka.stats.txt
  """
}

//CG Step: Perform core genome alignment using Roary
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

//CG Step: Infer ML tree from core genome alignment using IQ-TREE
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

//Run SNP calling pipeline if a reference sequence is provided
if (params.snp_reference != null & !params.snp_reference.isEmpty() | params.test_snp & params.test) {

    //SNP Step: Run CFSAN-SNP Pipeline
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

    //SNP Step: Run IQ-TREE on SNP alignment
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

    //SNP Step: Map cleaned reads to reference sequence using BWA
    process bwa {
      tag "$name"
      publishDir "${params.outdir}/mapping/sams", mode: 'copy',pattern:"*.sam"

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
  publishDir "${params.outdir}/mapping/bams", mode: 'copy',pattern:"*.sorted.bam*"
  publishDir "${params.outdir}/mapping/depth", mode: 'copy',pattern:"*.depth.tsv"
  publishDir "${params.outdir}/mapping/stats", mode: 'copy',pattern:"*.stats.txt"

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
  publishDir "${params.outdir}/mapping", mode: 'copy'

  input:
  file(depth) from assembly_depth_results.collect()
  file(mapped) from assembly_mapped_results.collect()

  output:
  file('coverage_stats.tsv') into assembly_mapping_tsv

  script:
  """
  #!/usr/bin/env python3
  import glob
  import os
  from numpy import median
  from numpy import average

  results = []
  files = glob.glob("*.assembly.depth.tsv*")
  for file in files:
    nums = []
    sid = os.path.basename(file).split('.')[0]
    with open(file,'r') as inFile:
      for line in inFile:
        nums.append(int(line.strip().split()[2]))
      med = int(median(nums))
      avg = int(average(nums))
      results.append(f"{sid}\\t{med}\\t{avg}\\n")
  with open('coverage_stats.tsv', 'w') as outFile:
    outFile.write("Sample\\tMedian Coverage (Mapped to Assembly)\\tMean Coverage (Mapped to Assembly)\\n")
    for result in results:
      outFile.write(result)
  """
}

if (params.snp_reference != null & !params.snp_reference.isEmpty() | params.test_snp & params.test) {

    //QC Step: Calculate coverage for reads mapped to reference
    process reference_mapping_stats {
      publishDir "${params.outdir}/mapping", mode: 'copy'

      input:
      file(depth) from reference_depth_results.collect()
      file(mapped) from reference_mapped_results.collect()

      output:
      file('mapping_results.tsv') into reference_mapping_tsv

      script:
      """
      #!/usr/bin/env python3
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
  // Set reference_mapping_tsv to non-tsv <- change this to an empty non-tsv file at some point
  reference_mapping_tsv = mapping_reference
//    reference_mapping_tsv = Channel.empty()
}

//QC Step: Merge QC results into one tsv
process merge_results {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  file(bbduk) from bbduk_tsv
  file(quast) from quast_tsv
  file(assembly) from assembly_mapping_tsv
  file(kraken) from kraken_tsv
  // make reference_mapping_tsv optional input
  //file(reference) from reference_mapping_tsv.ifEmpty{ 'empty' }
  file(reference) from reference_mapping_tsv

  output:
  file('dryad_report.csv')

  script:
  """
  #!/usr/bin/env python3

  import os
  import glob
  import pandas as pd
  from functools import reduce

  files = glob.glob('*.tsv')

  dfs = []

  for file in files:
      df = pd.read_csv(file, header=0, delimiter='\\t')
      dfs.append(df)

  merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],
                                              how='left'), dfs)

  merged = merged.rename(columns={'Contigs':'Contigs (#)'})

  merged.to_csv('dryad_report.csv', index=False, sep=',', encoding='utf-8')
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
//  file(c) from stats_multiqc.collect()
  file(e) from kraken_multiqc.collect()
  file(f) from logo
  file(config) from multiqc_config
  file(d) from prokka_multiqc.collect()
  file(r) from fastqc_multiqc.collect()

  output:
  file("*.html") into multiqc_output

  script:
  """
  multiqc -c ${config} .
  """
}
