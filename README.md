## Dryad is a bioinformatic pipeline to construct reference-free or SNP phylogenetic trees for examining prokarote relatedness in outbreaks.

Dryad uses multiprocessing and is able to construct a core genome tree and a SNP tree in parallel. Each process is set to use 4 cores. Running both tree construction pipelines would require 8 cores total.

#### Core Genome phylogenetic tree construction
The core genome tree pipeline takes a list of assembled genome locations and runs the following programs to construct a maximum likelihood tree.

Prokka v1.12 (https://github.com/tseemann/prokka)
Prokka is used to annotate the genomes.

Roary v3.6.0 (https://github.com/sanger-pathogens/Roary)
Roary takes the annotated genomes and constructs a core gene alignment.

RAxML v8.0.0 (https://sco.h-its.org/exelixis/web/software/raxml/index.html)
RAxML uses the core gene alignment and creates a maximum likelihood phylogenetic tree bootstraped 1000 times.

#### SNP phylogenetic tree construction
The SNP tree pipeline takes a list of paired end fastq file locations and a reference genome to construct a SNP based maximum likelihood tree.

The SNP pipeline uses Trimmomatic v0.35 (http://www.usadellab.org/cms/?page=trimmomatic) to trim the raw reads and Lyve-SET v2.0 (https://github.com/lskatz/lyve-SET) to construct the tree.

#### Dependencies
Python 3 (https://www.python.org)
Docker (https://www.docker.com)
python Docker library (https://github.com/docker/docker-py)
