#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/prs
========================================================================================
 lifebit-ai/prs Polygenic Risk Scores Pipeline to calculate genetic risk score for any given phenotype/trait (using PRSice)
 #### Homepage / Documentation
 https://github.com/lifebit-ai/prs
----------------------------------------------------------------------------------------
*/

// Initialise variable to store optional parameters
extra_flags = ""

// Base Dataset or Discovery Dataset
Channel
  .fromPath(params.base)
  .ifEmpty { exit 1, "GWAS association file not found: ${params.base}" }
  .set { base }
index = params.index ? 'T' : 'F'
if ( params.A1 ) { extra_flags += " --A1 ${params.A1}" }
if ( params.A2 ) { extra_flags += " --A2 ${params.A2}" }
if ( params.chr ) { extra_flags += " --chr ${params.chr}" }
if ( params.stat ) { extra_flags += " --stat ${params.stat}" }
if ( params.snp ) { extra_flags += " --snp ${params.snp}" }
if ( params.bp ) { extra_flags += " --bp ${params.bp}" }
if ( params.pvalue ) { extra_flags += " --pvalue ${params.pvalue}" }

// Target dataset
Channel
  .fromFilePairs("${params.target}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
  .ifEmpty { error "No plink files matching: ${params.target}.{bed,bim,fam}" }
  .set { plink_targets }
Channel
  .fromPath(params.pheno)
  .ifEmpty { exit 1, "Phenotype file not found: ${params.pheno}" }
  .set { pheno }
ignore_fid = params.ignore_fid ? 'T' : 'F'
if ( params.pheno_col ) { extra_flags += " --pheno-col ${params.pheno_col}"}

// Clumping
no_clump = params.no_clump ? 'T' : 'F'
if ( params.proxy ) { extra_flags += " --proxy ${params.proxy}" }
// LD
Channel
  .fromPath(params.ld)
  .ifEmpty { exit 1, "LD reference file not found: ${params.ld}" }
  .set { ld }
if ( params.ld_dose_thres ) { extra_flags += " --ld-dose-thres ${params.ld_dose_thres}" }
if ( params.ld_geno ) { extra_flags += " --ld-geno ${params.ld_geno}" }
if ( params.ld_info ) { extra_flags += " --ld-info ${params.ld_info}"}
if ( params.ld_maf ) { extra_flags += " --ld-maf ${params.ld_maf}"}

// Covariate
Channel
  .fromPath(params.cov)
  .ifEmpty { exit 1, "Covariate file not found: ${params.cov}" }
  .set { cov }
if ( params.cov_col ) { extra_flags += " --cov-col ${params.cov_col}"}

// Polygenic Risk Calculations
if ( params.no_regress ) { extra_flags += " --no-regress"}
if ( params.all_score ) { extra_flags += " --all-score"}
if ( params.perm ) { extra_flags += " --perm ${params.perm}"}
if ( params.print_snp ) { extra_flags += " --print-snp"}

// R Markdown report
Channel
  .fromPath(params.rmarkdown)
  .ifEmpty { exit 1, "R Markdown script not found: ${params.rmarkdown}" }
  .set { rmarkdown  }
Channel
  .fromPath(params.quantile_plot)
  .ifEmpty { exit 1, "TXT quantile plot file for Rmd not found: ${params.quantile_plot}" }
  .set { quantile_plot  }

/*--------------------------------------------------
  Polygenic Risk Calculations
---------------------------------------------------*/

process polygen_risk_calcs {
  tag "$name"
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file base from base
  set val(name), file(bed), file(bim), file(fam) from plink_targets
  file pheno from pheno
  file ld from ld
  file cov from cov

  output:
  file('*') into results
  file('*.png') into plots

  shell:
  pheno_flag = params.pheno.endsWith("no_pheno.txt") ? '' : "--pheno ${pheno}"
  ld_flag    = params.ld.endsWith("no_ld.txt")       ? '' : "--ld ${ld}"
  cov_flag   = params.cov.endsWith("no_cov.txt")     ? '' : "--cov ${cov}"
  '''
  PRSice.R \
    --prsice /usr/local/bin/PRSice_linux \
    --base !{base} \
    --index !{index} \
    --target !{name} \
    --ignore-fid !{ignore_fid} \
    --thread !{task.cpus} \
    --clump-kb !{params.clump_kb} \
    --clump-r2 !{params.clump_r2} \
    --clump-p !{params.clump_p} \
    --no-clump !{no_clump} \
    --missing !{params.missing} \
    --ld-hard-thres !{params.ld_hard_thres} \
    --model !{params.model} \
    --score !{params.score} \
    --quantile !{params.quantile} !{pheno_flag} !{ld_flag} !{cov_flag} !{extra_flags}

  # remove date from image names
  images=$(ls *.png)
  for image in $images; do
    date=$(echo $image | grep -Eo '_[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
    if ! [[ -z "${date// }" ]]; then
      mv "${image}" "${image/${date}/}"
    fi
  done
  '''
}

/*--------------------------------------------------
  Produce R Markdown report
---------------------------------------------------*/

process produce_report {
  publishDir params.outdir, mode: 'copy'

  input:
  file plots from plots
  file rmarkdown from rmarkdown
  file quantile_plot from quantile_plot

  output:
  file('*') into reports

  script:
  quantile_cmd = params.quantile ? "cat $quantile_plot >> $rmarkdown" : ''
  """
  # copy the rmarkdown into the pwd
  cp $rmarkdown tmp && mv tmp $rmarkdown

  $quantile_cmd

  R -e "rmarkdown::render('${rmarkdown}')"
  mkdir MultiQC && mv ${rmarkdown.baseName}.html MultiQC/multiqc_report.html
  """
}
