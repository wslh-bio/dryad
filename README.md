![Dryad](/assets/dryad_logo_250.png)

## Note: Dryad is currently undergoing large changes and may not be functional

Dryad is a [NextFlow](https://www.nextflow.io/) pipeline to construct reference free core-genome or SNP phylogenetic trees for examining prokaryote relatedness in outbreaks. Dryad will performs both a reference free core-genome analysis based off of the approach outlined by [Oakeson et. al](https://www.ncbi.nlm.nih.gov/pubmed/30158193) and/or a SNP analysis using the [CFSAN-SNP](https://snp-pipeline.readthedocs.io/en/latest/readme.html) pipeline.

### Table of Contents:
[Usage](#using-the-pipeline)  
[Workflow outline](#workflow-outline)  
[Core-genome](#core-Genome-phylogenetic-tree-construction)  
[SNP](#snp-phylogenetic-tree-construction)  
[Quality assessment](#quality-assessment)  
[Genome cluster report](#genome-cluster-report)  
[Output](#output-files)  
[Dependencies](#dependencies)  

### Using the pipeline
The pipeline is designed to start from raw Illumina short reads. All reads must be in the same directory. Then start the pipeline using `nextflow run k-florek/dryad`.

### Workflow outline

![Workflow](/assets/dryad_workflow_2.0.0.png)

Both pipelines begin with a quality trimming step to trim the reads of low quality bases at the end of the read using [Trimmomatic v0.39](http://www.usadellab.org/cms/?page=trimmomatic), the removal of PhiX contamination using [BBtools v38.76](https://jgi.doe.gov/data-and-tools/bbtools/), and the assessment of read quality using [FastQC v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). After processing, the reads are used by each pipeline as needed.  
*Note: Both pipelines can be run automatically in succession using the -cg and -s parameters simultaneously.*

### Core Genome phylogenetic tree construction
The core genome pipeline takes the trimmed and cleaned reads and infers a phylogenetic tree that can be used for inferring outbreak relatedness. This pipeline is based loosely off of the pipeline described here by [Oakeson et. al](https://www.ncbi.nlm.nih.gov/pubmed/30158193).

Species and MLST type are predicted from the assemblies generated during the core genome pipeline, and assembly quality is evaluated.

Additionally, the core genome pipeline can be run with `-ar` to predict antibiotic resistance genes.

The core genome pipeline uses the following applications and pipelines:

[Shovill v1.0.4](https://github.com/tseemann/shovill)
Shovill is a pipeline centered around SPAdes but alters some of the steps to get similar results in less time.

[Prokka v1.14.5](https://github.com/tseemann/prokka)
Prokka is a whole genome annotation tool that is used to annotate the coding regions of the assembly.

[Roary v3.12.0](https://github.com/sanger-pathogens/Roary)
Roary takes the annotated genomes and constructs a core gene alignment.

[IQ-Tree v1.6.7](http://www.iqtree.org/)
IQ-Tree uses the core gene alignment and creates a maximum likelihood phylogenetic tree bootstraped 1000 times.

[Mash v2.1](https://github.com/marbl/Mash)
Mash performs fast genome and metagenome distance estimation using MinHash.

[MLST v2.17.6](https://github.com/tseemann/mlst)
MLST scans contig files against PubMLST typing schemes.

[QUAST v5.0.2](http://bioinf.spbau.ru/quast)
QUAST evaluates genome assemblies.

[AMRFinderPlus v3.1.1](https://github.com/ncbi/amr)
AMRFinderPlus identifies acquired antimicrobial resistance genes.

### SNP phylogenetic tree construction
The SNP pipeline takes the trimmed and cleaned reads and infers a phylogenetic tree that can be used for inferring outbreak relatedness. The pipeline requires the path to the raw reads (mentioned above) and a reference genome in fasta file format.

The SNP pipeline uses the following applications and pipelines:

[CFSAN SNP Pipeline v2.0.2](https://github.com/CFSAN-Biostatistics/snp-pipeline)

[IQ-Tree v1.6.7](http://www.iqtree.org/)
IQ-Tree uses an alignment of the SNP sites to create a maximum likelihood phylogenetic tree bootstrapped 1000 times.

### Quality Assessment
The results of quality checks from each pipeline are summarized using [MultiQC v1.8](https://multiqc.info/)

### Genome cluster report
Dryad can generate an easily attributable analysis report. This uses RMarkdown and the results from the SNP and core genome pipelines to generate the genome cluster report. This option can be run using `--report`. The plotting defaults of the RMarkdown file (/report/report.Rmd) can be modified as necessary and rebuilt using `dryad_report`.

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
**\*.gff** - Annotations predicted by Prokka  
**\*.prokka.stats.txt** - Prokka log files  
**core_gene_alignment.aln** - Core-genome alignment  
**\*.contigs.fa** - Shovill assemblies  
**core_genome_statistics.txt** - Text file with the number of genes in the core and accessory genomes  
**core_genome.tree** - ML tree inferred from core-genome alignment  
**dryad_report.csv** - Summary table of each step in Spriggan  
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
[Kelsey Florek](https://github.com/k-florek), WSLH Bioinformatics Scientist  
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatics Scientist
