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
  file('*') into plots

  script:
  """
  PRSice_v1.2.R -q --args \
  plink /usr/local/bin/plink \
  base $assoc  \
  target $name \
  quantiles $quantiles
  """
}
