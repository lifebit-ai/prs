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


log.info """\
L I F E B I T - A I / P R S  P I P E L I N E
==========================================================
saige_base                : ${params.saige_base}
gwas_catalogue_base       : ${params.gwas_catalogue_base}
target_plink_files_dir    : ${params.target_plink_dir}
target_pheno              : ${params.target_pheno}
outdir                    : ${params.outdir}             
"""


/*-----------------------------------------------------
  Setting up base dataset: SAIGE or GWAS catalogue data
------------------------------------------------------*/

// Initialise variable to store optional parameters
extra_flags = ""

// Base Dataset or Discovery Dataset

if (params.saige_base) {
  saige_base_ch = Channel
      .fromPath(params.saige_base, checkIfExists: true)
      .ifEmpty { exit 1, "SAIGE summary stats (base cohort) not found: ${params.saige}" }
  
  extra_flags += " --A1 Allele1"
  extra_flags += " --A2 Allele2"
  extra_flags += " --chr CHR"
  extra_flags += " --stat BETA"
  extra_flags += " --snp SNPID"
  extra_flags += " --bp POS"
  extra_flags += " --pvalue p.value"
  
} else if (params.gwas_catalogue_base){
  // gwas_catatologue_base_ch = Channel
  //    .fromPath(params.gwas_catalogue_base, checkIfExists: true)
  //    .ifEmpty { exit 1, "GWAS summary stats (base cohort) not found: ${params.base}" }
  
  // THESE FLAGS WILL BE PREDICTABLE ONCE I HAVE SORTED GWAS CATALOGUE
  // extra_flags += " --A1 Allele1"
  // extra_flags += " --A2 Allele2"
  // extra_flags += " --chr CHR"
  // extra_flags += " --stat BETA"
  // extra_flags += " --snp SNPID"
  // extra_flags += " --bp POS"
  // extra_flags += " --pvalue p.value"

}
// saige_base_ch.view()


/*-------------------------------------------------------------------------------------
  Setting up base dataset: transforming SAIGE output from lifebit-ai/gel-gwas
---------------------------------------------------------------------------------------*/

process transform_saige_base {
    publishDir "${params.outdir}", mode: "copy"

    input:
    file saige_base from saige_base_ch

    output:
    file("*") into transformed_saige_base_ch
    
    script:
    """
    transform_base_saige.R ${saige_base}
    """

}
transformed_saige_base_ch.view()


/*-------------------------------------------------------------------------------------
  Setting up target dataset: PLINK files and pheno file output from lifebit-ai/gel-gwas
---------------------------------------------------------------------------------------*/

if (params.target_plink_dir) {
    Channel
    .fromPath("${params.target_plink_dir}/*.{bed,bim,fam}")
    .ifEmpty { error "No target plink files found in : ${params.target_plink_dir}" }
    .set { target_plink_dir_ch}
}
// target_plink_files_ch.view()

Channel
  .fromPath(params.target_pheno, checkIfExists: true)
  .ifEmpty { exit 1, "Phenotype file not found: ${params.target_pheno}" }
  .set { target_pheno_ch }
// target_pheno_ch.view()

process transform_target_pheno {
    publishDir "${params.outdir}", mode: "copy"

    input:
    file pheno from target_pheno_ch

    output:
    file("*") into transformed_target_pheno_ch
    
    script:
    """
    transform_target_pheno.R ${pheno}
    """

}
//transformed_target_pheno_ch.view()


