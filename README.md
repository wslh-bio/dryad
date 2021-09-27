![Dryad](/assets/dryad_logo_250.png)

![Dryad](https://github.com/wslh-bio/dryad/actions/workflows/dryad_build.yml/badge.svg)
![GPL-3.0](https://img.shields.io/github/license/wslh-bio/dryad)
![GitHub Release](https://img.shields.io/github/release/wslh-bio/dryad)

Dryad is a [NextFlow](https://www.nextflow.io/) pipeline to construct reference free core-genome or SNP phylogenetic trees for examining prokaryote relatedness in outbreaks. Dryad will performs both a reference free core-genome analysis based off of the approach outlined by [Oakeson et. al](https://www.ncbi.nlm.nih.gov/pubmed/30158193) and/or a SNP analysis using the [CFSAN-SNP](https://snp-pipeline.readthedocs.io/en/latest/readme.html) pipeline.

### Table of Contents:
[Usage](#using-the-pipeline)  
[Workflow outline](#workflow-outline)  
[Core-genome](#core-genome-alignment-and-phylogenetic-tree-construction)  
[SNP](#snp-distances-calculation-and-phylogenetic-tree-construction)  
[Quality assessment](#quality-assessment)  
[Output](#output-files)  

### Using the pipeline
The pipeline is designed to start from raw Illumina short reads. All reads must be in the same directory. Then start the pipeline using:  
```
nextflow wslh-bio/dryad -r <version> --reads [path-to-reads]
```  
to run the SNP pipeline include the `--snp_reference` parameter:  
```
nextflow wslh-bio/dryad -r <version> --reads [path-to-reads] --snp_reference [path-to-reference-fasta]
```  

You can also test the pipeline with example data using `--test`, note this requires NextFlow version `21.07.0-edge` or greater:
```
nextflow dryad.nf --test
```

### Workflow outline

![Workflow](/assets/dryad_workflow_3.0.png)

### Read trimming and cleaning
Read trimming and cleaning is performed using [BBtools v38.76](https://jgi.doe.gov/data-and-tools/bbtools/) to trim reads of low quality bases and remove PhiX contamination. After processing, the reads are used by each pipeline as needed.  
*Note: Both pipelines can be run automatically in parallel by supplying the snp_reference parameter.*

### Core genome alignment and phylogenetic tree construction
The core genome pipeline takes the trimmed and cleaned reads and infers a phylogenetic tree that can be used for inferring outbreak relatedness. This pipeline is based loosely off of the pipeline described here by [Oakeson et. al](https://www.ncbi.nlm.nih.gov/pubmed/30158193).

Species is predicted from the cleaned reads, and assembly quality is evaluated.

The core genome pipeline uses the following applications and pipelines:

[Shovill v1.0.4](https://github.com/tseemann/shovill)
Shovill is a genome assembly pipeline that uses SPAdes, but it alters certain steps to get similar results in less time.

[Prokka v1.14.5](https://github.com/tseemann/prokka)
Prokka is a whole genome annotation tool used to annotate the coding regions of the assembly.

[Roary v3.12.0](https://github.com/sanger-pathogens/Roary)
Roary takes the annotated genomes and constructs a core gene alignment.

[IQ-Tree v1.6.7](http://www.iqtree.org/)
IQ-Tree uses the core gene alignment and infers a maximum likelihood phylogenetic tree bootstrapped 1000 times.

[QUAST v5.0.2](http://bioinf.spbau.ru/quast)
QUAST evaluates genome assembly quality.

[Kraken2 v2.0.8](https://ccb.jhu.edu/software/kraken2/)
Kraken uses the cleaned reads to predict species and detect contamination.

### SNP distances calculation and phylogenetic tree construction
The SNP pipeline takes the trimmed and cleaned reads and infers a phylogenetic tree that can be used for inferring outbreak relatedness. The pipeline requires the path to the raw reads (mentioned above) and a reference genome in fasta file format.

The SNP pipeline uses the following applications and pipelines:

[CFSAN SNP Pipeline v2.0.2](https://github.com/CFSAN-Biostatistics/snp-pipeline)

[IQ-Tree v1.6.7](http://www.iqtree.org/)
IQ-Tree uses an alignment of the SNP sites to create a maximum likelihood phylogenetic tree bootstrapped 1000 times.

#### Quality Assessment

[FastQC v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
FastQC is used assess the quality of the raw and cleaned reads.

[QUAST v5.0.2](http://bioinf.spbau.ru/quast)
QUAST assesses the quality of the genome assemblies.

[Samtools v1.10](http://www.htslib.org/)
calculates the number and depth of cleaned reads mapped to their assemblies and the reference genome. [BWA v0.7.17-r1188](http://bio-bwa.sourceforge.net/) is used for mapping reads.

[MultiQC v1.8](https://multiqc.info/)
summarizes the results of FastQC, Prokka, Samtools Stats and Kraken.
### Output files
```
dryad_results
├── annotated
│   ├── *.gff
│   └── *.prokka.stats.txt
├── assembled
│   └── *.contigs.fa
├── core_gene_alignment.aln
├── core_genome_statistics.txt
├── core_genome.tree
├── dryad_report.csv
├── assembled
├── fastqc
│   ├── fastqc_summary.txt
│   ├── *.html
│   └── zips
│       └── *.zip
├── kraken
│   ├── kraken_results.tsv
│   └── *kraken2_report.txt
├── mapping
│   ├── bams
│   │   ├── *.assembly.bam
│   │   ├── *.assembly.bai
│   │   ├── *.reference.bam
│   │   └── *.reference.bai
│   ├── coverage_stats.tsv
│   ├── depth
│   │   ├── *.assembly.depth.tsv
│   │   └── *.reference.depth.tsv
│   ├── mapping_stats.tsv
│   ├── sams
│   │   ├── *.assembly.sam
│   │   └── *.reference.sam
│   └── stats
│       ├── *.assembly.stats.txt
│       └── *.reference.stats.txt
├── multiqc_report.html
├── quast
│   ├── quast_results.tsv
│   └── *.quast.tsv
├── snp_distance_matrix.tsv
├── snpma.fasta
├── snp.tree
└── trimming
    ├── bbduk_results.tsv
    └── *.trim.txt
```
**\*.gff** - Gene annotations predicted by Prokka  
**\*.prokka.stats.txt** - Prokka log files  
**core_gene_alignment.aln** - Core-genome alignment  
**\*.contigs.fa** - Shovill assemblies  
**core_genome_statistics.txt** - Text file with the number of genes in the core and accessory genomes  
**core_genome.tree** - ML tree inferred from core-genome alignment  
**dryad_report.csv** - Summary table of each step in dryad  
**fastqc_summary.txt** - Summary table of FastQC results  
**\*.html** - HTML files of FastQC results  
**\*.zip** - FastQC results, compressed  
**kraken_results.tsv** - Summary table of Kraken results  
**\*kraken2_report.txt** - Report of Kraken results for each sample  
**\*.assembly.bam** - Alignments to an assembly BAM format  
**\*.assembly.bai** - Index file of alignments to an assembly  
**\*.reference.bam**  - Alignments to the reference sequence BAM format  
**\*.reference.bai** - Index file of alignments to the reference sequence  
**coverage_stats.tsv** - Summary table of mean and median coverage calculated with Samtools depth  
**\*.assembly.depth.tsv** - Raw Samtools depth output for reads mapped to an assembly  
**\*.reference.depth.tsv** - Raw Samtools depth output for reads mapped to the reference sequence  
**mapping_stats.tsv** - Summary table of statistics   
**\*.assembly.sam** - Alignments to an assembly SAM format  
**\*.reference.sam** - Alignments to the reference sequence SAM format  
**\*.assembly.stats.txt** - Output of Samtools stats for reads mapped to an assembly  
**\*.reference.stats.txt** - Output of Samtools stats for reads mapped to the reference sequence  
**multiqc_report.html** - HTML report generated by MultiQC  
**quast_results.tsv** - Summary table of QUAST results  
**\*.quast.tsv** - QUAST results for each sample  
**snp_distance_matrix.tsv** - SNP distance matrix  
**snpma.fasta** - SNP alignment  
**snp.tree** - SNP tree  
**bbduk_results.tsv** - Summary table of trimming with BBduk  
**\*.trim.txt** - Trimming results from BBduk each sample  

### Authors
[Kelsey Florek](https://github.com/k-florek), WSLH Senior Genomics and Data Scientist  
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatics Scientist
