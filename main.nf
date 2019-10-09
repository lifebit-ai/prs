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
  .fromPath(params.assoc)
  .ifEmpty { exit 1, "GWAS association file not found: ${params.assoc}" }
  .set { assoc }

// Target dataset
Channel
  .fromFilePairs("${params.target}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
  .ifEmpty { error "No plink files matching: ${params.target}.{bed,bim,fam}" }
  .set { plink_targets }

// Polygenic Risk Calculations
quantiles = params.quantiles ? 'T' : 'F'

// R Markdown report
Channel
  .fromPath(params.rmarkdown)
  .ifEmpty { exit 1, "R Markdown script not found: ${params.rmarkdown}" }
  .set { rmarkdown  }

/*--------------------------------------------------
  Polygenic Risk Calculations
---------------------------------------------------*/

process polygen_risk_calcs {
  tag "$name"
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file assoc from assoc
  set val(name), file(bed), file(bim), file(fam) from plink_targets

  output:
  file('*') into results
  file('*.png') into plots

  shell:
  '''
  PRSice_v1.2.R -q --args \
  plink /usr/local/bin/plink \
  base !{assoc}  \
  target !{name} \
  quantiles !{quantiles}

  # remove date from image names
  images=$(ls *.png)
  for image in $images; do
    date=$(echo $image | grep -Eo '[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
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
  file(plots) from plots
  file(rmarkdown) from rmarkdown

  output:
  file('*') into reports

  script:
  if (params.quantiles) {
    quantile_plot = """
                    Row
                    -------------------------------------
                        
                    ### Quantile Plot

                    ![PRSice_QUANTILE_PLOT](PRSice_QUANTILE_PLOT.png)

                    > Generated using PRSice
                    """
  } else { quantile_plot = '' }
  """
  # copy the rmarkdown into the pwd
  cp $rmarkdown tmp && mv tmp $rmarkdown

  echo $quantile_plot >> $rmarkdown

  R -e "rmarkdown::render('${rmarkdown}')"
  mkdir MultiQC && mv ${rmarkdown.baseName}.html MultiQC/multiqc_report.html
  """
}
