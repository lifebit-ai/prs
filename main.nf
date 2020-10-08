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
gwas_catalogue_base       : ${params.gwas_cat_study_id}
pheno_metadata            : ${params.pheno_metadata}
target_plink_files_dir    : ${params.target_plink_dir}
target_pheno              : ${params.target_pheno}
outdir                    : ${params.outdir}             
"""



/*-----------------------------------------------------
  Setting up base dataset: SAIGE or GWAS catalogue data
------------------------------------------------------*/

// Base Dataset or Discovery Dataset

if (params.saige_base) {
  saige_base_ch = Channel
      .fromPath(params.saige_base, checkIfExists: true)
      .ifEmpty { exit 1, "SAIGE summary stats (base cohort) not found: ${params.saige}" }
} else if (params.gwas_cat_study_id){
  gwas_catalogue_ftp_ch = Channel
      .fromPath(params.gwas_catalogue_ftp, checkIfExists: true)
      .ifEmpty { exit 1, "GWAS catalogue ftp locations not found: ${params.gwas_catalogue_ftp}" }
      .splitCsv(header: true)
      .map { row -> tuple(row.study_accession, row.ftp_link_harmonised_summary_stats) }
      .filter{ it[0] == params.gwas_cat_study_id}
      .ifEmpty { exit 1, "The GWAS study accession number you provided does not come as a harmonized dataset that can be used as a base cohort "}
      .flatten()
      .last()
} else {
  exit 1, "No SAIGE base or GWAS catalogue study was provided as base for this PRS!"
}


/*-------------------------
  Transforming SAIGE input 
---------------------------*/

if (params.saige_base) {
  process transform_saige_base {
    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"

    input:
    file saige_base from saige_base_ch

    output:
    file("*") into transformed_base_ch
    
    script:
    """
    transform_base_saige.R ${saige_base}
    """
    } 
}



/*---------------------------------
  Transforming GWAS catalogue input 
-----------------------------------*/

process download_gwas_catalogue {
  label 'high_memory'
  publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"
  
  input:
  val(ftp_link) from gwas_catalogue_ftp_ch

  output:
  file("*.h.tsv.gz") into downloaded_gwas_catalogue_ch

  script:
  """
  wget ${ftp_link}
  """
}


process transform_gwas_catalogue_base {
  label 'high_memory'
  publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"
  
  input:
  file gwas_catalogue_base from downloaded_gwas_catalogue_ch

  output:
  file("*") into transformed_base_ch

  script:
  """
  cp /opt/bin/* .

  transform_base_gwas_catalogue.R ${gwas_catalogue_base}
  """
}

/*-------------------------------------------------------------------------------------
  Obtaining phenotype metadata - necessary for determining which covariates to plot
--------------------------------------------------------------------------------------*/

/* pheno_metadata_ch = Channel
      .fromPath(params.pheno_metadata, checkIfExists: true)
      .ifEmpty { exit 1, "Phenotype metadata file not found: ${params.pheno_metadata}" }

 */

/*-------------------------------------------------------------------------------------
  Setting up target dataset: PLINK files and pheno file output from lifebit-ai/gel-gwas
---------------------------------------------------------------------------------------*/

/* if (params.target_plink_dir) {
    Channel
    .fromPath("${params.target_plink_dir}/*.{bed,bim,fam}")
    .ifEmpty { error "No target plink files found in : ${params.target_plink_dir}" }
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple()
    .set { target_plink_dir_ch }
}

Channel
  .fromPath(params.target_pheno, checkIfExists: true)
  .ifEmpty { exit 1, "Phenotype file not found: ${params.target_pheno}" }
  .set { target_pheno_ch }

/*---------------------------------
  Transforming target pheno input 
-----------------------------------*/

/* process transform_target_pheno {
    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"

    input:
    file pheno from target_pheno_ch

    output:
    file("*") into transformed_target_pheno_ch 
    file("target.pheno") into transformed_target_pheno_for_plots_ch
    file("target.cov") into transformed_target_cov_for_plots_ch
    
    script:
    """
    transform_target_pheno.R ${pheno}
    """

} */

/*----------------------------
  Setting up other parameters
------------------------------*/

// Clumping

/* no_clump = params.no_clump ? 'T' : 'F'
if ( params.proxy ) { extra_flags += " --proxy ${params.proxy}" } */

// LD

/* Channel
  .fromPath(params.ld)
  .ifEmpty { exit 1, "LD reference file not found: ${params.ld}" }
  .set { ld } */

/* if ( params.ld_dose_thres ) { extra_flags += " --ld-dose-thres ${params.ld_dose_thres}" }
if ( params.ld_geno ) { extra_flags += " --ld-geno ${params.ld_geno}" }
if ( params.ld_info ) { extra_flags += " --ld-info ${params.ld_info}"}
if ( params.ld_maf ) { extra_flags += " --ld-maf ${params.ld_maf}"} */

// Polygenic Risk Calculations

/* if ( params.no_regress ) { extra_flags += " --no-regress"}
if ( params.all_score ) { extra_flags += " --all-score"}
if ( params.perm ) { extra_flags += " --perm ${params.perm}"}
if ( params.print_snp ) { extra_flags += " --print-snp"} */

// R Markdown report

/* Channel
  .fromPath(params.rmarkdown)
  .ifEmpty { exit 1, "R Markdown script not found: ${params.rmarkdown}" }
  .set { rmarkdown  }

Channel
  .fromPath(params.quantile_plot)
  .ifEmpty { exit 1, "TXT quantile plot file for Rmd not found: ${params.quantile_plot}" }
  .set { quantile_plot  } */


/*--------------------------------------------------
  Polygenic Risk Calculations
---------------------------------------------------*/

/* process polygen_risk_calcs {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file base from transformed_base_ch
  tuple val(name), file("*") from target_plink_dir_ch
  tuple file(pheno), file(cov) from transformed_target_pheno_ch

  output:
  file("*") into results
  file('PRSice.best') into results_for_plots
  file('*.png') into plots_p1

  shell:
  '''
  PRSice.R \\
    --prsice /usr/local/bin/PRSice_linux \\
    --base !{base} \\
    --snp SNPID \\
    --chr CHR \\
    --bp POS \\
    --A1 Allele1 \\
    --A2 Allele2 \\
    --stat BETA \\
    --pvalue p.value \\
    --beta \\
    --target !{name}_chr#_filtered \\
    --binary-target !{params.binary_trait} \\
    --pheno !{pheno} \\
    --cov !{cov} \\
    --thread !{task.cpus} \\
    --clump-kb !{params.clump_kb} \\
    --clump-r2 !{params.clump_r2} \\
    --clump-p !{params.clump_p} \\
    --no-clump !{no_clump} \\
    --missing !{params.missing} \\
    --ld-hard-thres !{params.ld_hard_thres} \\
    --model !{params.model} \\
    --score !{params.score} \\
    --quantile !{params.quantile} \\

  # Remove date from image names (only for images produced by PRSice)
  images=$(ls *.png)
  for image in $images; do
    date=$(echo $image | grep -Eo '_[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
    if ! [[ -z "${date// }" ]]; then
      mv "${image}" "${image/${date}/}"
    fi
  done
   '''
} */

/*--------------------------------------------------
  Additional visualizations
---------------------------------------------------*/

/* process additional_plots {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file pheno from transformed_target_pheno_for_plots_ch
  file cov from transformed_target_cov_for_plots_ch
  file prs from results_for_plots
  file metadata from pheno_metadata_ch

  output:
  file('*.png') into plots_p2

  script:
  """
  cp /opt/bin/* .

  plot_prs_vs_cov.R ${cov} ${prs} ${metadata}
  plot_prs_vs_density.R ${pheno} ${prs}
  """

} */

/*--------------------------------------------------
  Produce R Markdown report                          
---------------------------------------------------*/

// Concatenate plot channels
/* all_plots_ch = plots_p1.concat(plots_p2).flatten().toList()

process produce_report {
  publishDir params.outdir, mode: 'copy'

  input:
  file plots from all_plots_ch
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
 */