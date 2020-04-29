#!/usr/bin/env nextflow

//Description: Workflow for regenerating report
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

params.snp_matrix = ""
params.cg_tree = ""
params.ar_tsv = ""
params.rmd = ""

Channel.fromPath(params.snp_matrix).set{ snp_mat }
Channel.fromPath(params.cg_tree).set{ cgtree }
Channel.fromPath(params.ar_tsv).set{ ar_tsv }
Channel.fromPath(params.report).set{ report }

process render{
  publishDir "${params.outdir}/report", mode: 'copy', pattern: "*.[pdf,Rmd]"

  input:
  file snp from snp_mat
  file tree from cgtree
  file ar from ar_tsv
  file rmd from report

  output:
  file "cluster_report.pdf"
  file "report_template.Rmd"

  shell:
  """
  Rscript /reports/render_dryad.R ${snp} ${tree} ${ar} ${rmd}
  mv report.pdf cluster_report.pdf
  """
}
