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

// Base Dataset or Discovery Dataset
Channel
  .fromPath(params.base)
  .ifEmpty { exit 1, "GWAS association file not found: ${params.base}" }
  .set { base }

// Target dataset
Channel
  .fromFilePairs("${params.target}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
  .ifEmpty { error "No plink files matching: ${params.target}.{bed,bim,fam}" }
  .set { plink_targets }

// Clumping
no_clump = params.no_clump ? 'T' : 'F'
proxy = params.proxy ? '--proxy !{params.proxy}' : ''

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

  output:
  file('*') into results
  file('*.png') into plots

  shell:
  '''
  PRSice.R \
    --prsice /usr/local/bin/PRSice_linux \
    --base !{base} \
    --target !{name} \
    --thread !{task.cpus} \
    --clump-kb !{params.clump_kb} \
    --clump-r2 !{params.clump_r2} \
    --clump-p !{params.clump_p} \
    --no-clump !{no_clump} \
    --missing !{params.missing} \
    --quantile !{params.quantile} !{proxy}

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
