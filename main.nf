#!/usr/bin/env nextflow
/*
===============================================================================
                         lifebit-ai/prs
===============================================================================
 Polygenic Risk Scores Pipeline to calculate genetic risk score for any  given 
 phenotype/trait (using PRSice)

 Homepage / Documentation: https://github.com/lifebit-ai/prs
-------------------------------------------------------------------------------
*/



log.info """\
L I F E B I T - A I / P R S  P I P E L I N E
==========================================================
saige_base                : ${params.saige_base}
gwas_cat_study_id         : ${params.gwas_cat_study_id}
pheno_metadata            : ${params.pheno_metadata}
target_plink_files_dir    : ${params.target_plink_dir}
target_pheno              : ${params.target_pheno}
binary_trait              : ${params.binary_trait}
outdir                    : ${params.outdir}             
"""



/*-------------------------------------------------------------------------------------
  Setting up target dataset: PLINK files and pheno file output from lifebit-ai/gel-gwas
---------------------------------------------------------------------------------------*/

if (params.target_plink_dir) {
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



/*-------------------------------------------------------------
  Setting up target pheno file: output from lifebit-ai/gel-gwas
---------------------------------------------------------------*/

Channel
  .fromPath(params.target_pheno, checkIfExists: true)
  .ifEmpty { exit 1, "Phenotype file not found: ${params.target_pheno}" }
  .set { target_pheno_ch }



/*---------------------------------
  Transforming target pheno input 
-----------------------------------*/

process transform_target_pheno {
    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"

    input:
    file pheno from target_pheno_ch

    output:
    tuple file("target.pheno"), file("target.cov") into (transformed_target_pheno_ch, transformed_target_pheno_for_plots_ch)

    script:
    """
    transform_target_pheno.R --input_pheno ${pheno}
    """

}



/*----------------------------------------------------------------------------------
  Setting up phenotype metadata: necessary for determining which covariates to plot
------------------------------------------------------------------------------------*/

pheno_metadata_ch = Channel
      .fromPath(params.pheno_metadata, checkIfExists: true)
      .ifEmpty { exit 1, "Phenotype metadata file not found: ${params.pheno_metadata}" }



/*-----------------------------------------------------
  Setting up base dataset: SAIGE or GWAS catalogue data
-------------------------------------------------------*/

// Base Dataset or Discovery Dataset

if (params.saige_base && params.gwas_cat_study_id) {
  exit 1, "You have provided both SAIGE summary statistics and GWAS catalogue summary statistics: only one base cohort can be used. \
  \nPlease use only one of these inputs (i.e. only --saige_base or only --gwas_cat_study_id)."
}

if (!params.saige_base && !params.gwas_cat_study_id) {
  exit 1, "No SAIGE summary statistics or GWAS catalogue summary statistics were provided as base for this PRS!  \
  \nPlease provide either a SAIGE output file (.csv) or a GWAS catalogue study accession id (for example GCST004420)."
}

if (params.saige_base) {
  saige_base_ch = Channel
      .fromPath(params.saige_base, checkIfExists: true)
      .ifEmpty { exit 1, "SAIGE summary stats (base cohort) not found: ${params.saige}" }
}

if (params.gwas_cat_study_id){
  gwas_catalogue_ftp_ch = Channel
      .fromPath(params.gwas_catalogue_ftp, checkIfExists: true)
      .ifEmpty { exit 1, "GWAS catalogue ftp locations not found: ${params.gwas_catalogue_ftp}" }
      .splitCsv(header: true)
      .map { row -> tuple(row.study_accession, row.ftp_link_harmonised_summary_stats) }
      .filter{ it[0] == params.gwas_cat_study_id}
      .ifEmpty { exit 1, "The GWAS study accession number you provided does not come as a harmonized dataset that can be used as a base cohort "}
      .flatten()
      .last()
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
    file("base.data") into transformed_base_ch
    
    script:
    """
    transform_base_saige.R --input_saige ${saige_base}
    """
    } 
}



/*---------------------------------
  Transforming GWAS catalogue input 
-----------------------------------*/

if (params.gwas_cat_study_id) {
  process download_gwas_catalogue {
    label "high_memory"
    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"
    
    input:
    val(ftp_link) from gwas_catalogue_ftp_ch
    
    output:
    file("*") into downloaded_gwas_catalogue_ch
    
    script:
    """
    wget ${ftp_link}
    """
  }

  process transform_gwas_catalogue_base {
    label "high_memory"
    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"
    
    input:
    file gwas_catalogue_base from downloaded_gwas_catalogue_ch
    
    output:
    file("base.data") into transformed_base_ch
    
    script:
    """
    transform_base_gwas_catalogue.R --input_gwas_cat ${gwas_catalogue_base}
    """
    }
}



/*----------------------------
  Setting up other parameters
------------------------------*/

// Initialise variable to store optional parameters
extra_flags = ""

// 1 - Base file options

if ( params.base_info ) { extra_flags += " --base-info ${params.base_info}" }
if ( params.base_maf ) { extra_flags += " --base-maf ${params.base_maf}" }
if ( params.no_default ) { extra_flags += " --no-default ${params.no_default}" }

// 2 - Target file options (other than the hardcoded ones in main.nf)

if ( params.geno ) { extra_flags += " --geno ${params.geno}" }
if ( params.info ) { extra_flags += " --info ${params.info}" }
if ( params.keep ) { extra_flags += " --keep ${params.keep}" }
if ( params.maf ) { extra_flags += " --maf ${params.maf}" }
if ( params.nonfounders ) { extra_flags += " --nonfounders ${params.nonfounders}" }
if ( params.ignore_fid ) { extra_flags += " --ignore-fid ${params.ignore_fid}" }
if ( params.prevalence ) { extra_flags += " --prevalence ${params.prevalence}" }
if ( params.remove ) { extra_flags += " --remove ${params.remove}" }

// 3 - Dosage commands not available (pipeline is not currently developed to handle dosages)

// 4 - Clumping

if ( params.clump_kb ) { extra_flags += " --clump-kb ${params.clump_kb}" }
if ( params.clump_r2 ) { extra_flags += " --clump-r2 ${params.clump_r2}" }
if ( params.clump_p ) { extra_flags += " --clump-p ${params.clump_p}" }
if ( params.ld ) { extra_flags += " --ld ${params.ld}" }
if ( params.ld_dose_thres ) { extra_flags += " --ld-dose-thres ${params.ld_dose_thres}" }
if ( params.ld_geno ) { extra_flags += " --ld-geno ${params.ld_geno}" }
if ( params.ld_info ) { extra_flags += " --ld-info ${params.ld_info}" }
if ( params.ld_hard_thres ) { extra_flags += " --ld-hard-thres ${params.ld_hard_thres}" }
if ( params.ld_keep ) { extra_flags += " --ld-keep ${params.ld_keep}" }
if ( params.ld_list ) { extra_flags += " --ld-list ${params.ld_list}" }
if ( params.ld_maf ) { extra_flags += " --ld-maf ${params.ld_maf}" }
if ( params.ld_remove ) { extra_flags += " --ld-remove ${params.ld_remove}" }
if ( params.ld_type ) { extra_flags += " --ld-type ${params.ld_type}" }
if ( params.no_clump ) { extra_flags += " --no-clump ${params.no_clump}" }
if ( params.proxy ) { extra_flags += " --proxy ${params.proxy}" }

// 5 - Covariate commands not available as pipeline handles the covariates directly

// 6 - P-value thresholding

if ( params.bar_levels ) { extra_flags += " --bar-levels ${params.bar_levels}" }
if ( params.fastscore ) { extra_flags += " --fastscore ${params.fastscore}" }
if ( params.no_full ) { extra_flags += " --no-full ${params.no_full}" }
if ( params.interval ) { extra_flags += " --interval ${params.interval}" }
if ( params.lower ) { extra_flags += " --lower ${params.lower}" }
if ( params.model ) { extra_flags += " --model ${params.model}" }
if ( params.missing ) { extra_flags += " --missing ${params.missing}" }
if ( params.no_regress ) { extra_flags += " --no-regress ${params.no_regress}" }
if ( params.score ) { extra_flags += " --score ${params.score}" }
if ( params.upper ) { extra_flags += " --upper ${params.upper}" }

// 7 - R specific commands not available (pipeline is using a Docker image that handles such details)

// 8 - Plotting

if ( params.bar_col_high ) { extra_flags += " --bar-col-high ${params.bar_col_high}" }
if ( params.bar_col_low ) { extra_flags += " --bar-col-low ${params.bar_col_low}" }
if ( params.bar_col_p ) { extra_flags += " --bar-col-p ${params.bar_col_p}" }
if ( params.bar_palatte ) { extra_flags += " --bar-palatte ${params.bar_palatte}" }
if ( params.multi_plot ) { extra_flags += " --multi-plot ${params.multi_plot}" }
if ( params.plot ) { extra_flags += " --plot ${params.plot}" }
if ( params.plot_set ) { extra_flags += " --plot-set ${params.plot_set}" }
if ( params.scatter_r2 ) { extra_flags += " --scatter-r2 ${params.scatter_r2}" }

// 8 - Miscellaneous

if ( params.all_score ) { extra_flags += " --all-score ${params.all_score}" }
if ( params.exclude ) { extra_flags += " --exclude ${params.exclude}" }
if ( params.extract ) { extra_flags += " --extract ${params.extract}" }
if ( params.ignore_fid ) { extra_flags += " --ignore-fid ${params.ignore_fid}" }
if ( params.keep_ambig ) { extra_flags += " --keep-ambig ${params.keep_ambig}" }
if ( params.logit_perm ) { extra_flags += " --logit-perm ${params.logit_perm}" }
if ( params.memory ) { extra_flags += " --memory ${params.memory}" }
if ( params.non_cumulate ) { extra_flags += " --non-cumulate ${params.non_cumulate}" }
if ( params.out ) { extra_flags += " --out ${params.out}" }
if ( params.perm ) { extra_flags += " --perm ${params.perm}" }
if ( params.print_snp ) { extra_flags += " --print-snp ${params.print_snp}" }
if ( params.seed ) { extra_flags += " --seed ${params.seed}" }
if ( params.thread ) { extra_flags += " --thread ${params.thread}" }
if ( params.ultra ) { extra_flags += " --ultra ${params.ultra}" }
if ( params.x_range ) { extra_flags += " --x-range ${params.x_range}" }



/*----------------------------
  Polygenic Risk Calculations
------------------------------*/

process polygen_risk_calcs {
  publishDir "${params.outdir}", mode: "copy"

  input:
  file base from transformed_base_ch
  tuple val(name), file("*") from target_plink_dir_ch
  tuple file(pheno), file(cov) from transformed_target_pheno_ch

  output:
  file("*") into all_results_ch
  file("PRSice.best") into best_PRS_ch
  file("PRSice*") into results_for_report_ch

  shell:
  quantile_flag = params.quantile =~ false ? '' : "--quantile ${params.quantile}"
  quant_break_flag = params.quant_break =~ false ? '' : "--quant-break ${params.quant_break}"
  quant_ref_flag = params.quant_ref =~ false ? '' : "--quant-ref ${params.quant_ref}"
  quant_extract_flag = params.quant_extract =~ false ? '' : "--quant-extract ${params.quant_extract}"
  '''
  PRSice.R \\
    --prsice /usr/local/bin/PRSice_linux \\
    --base !{base} \\
    --snp SNPID \\
    --chr CHR \\
    --bp POS \\
    --A1 Allele1 \\
    --A2 Allele2 \\
    --pvalue p.value \\
    --stat BETA \\
    --beta \\
    --target !{name}_chr#_filtered \\
    --binary-target !{params.binary_trait} \\
    --pheno !{pheno} \\
    --cov !{cov} !{extra_flags} !{quantile_flag} !{quant_break_flag} !{quant_ref_flag} !{quant_extract_flag} \\

  # Remove date from image names (only for images produced by PRSice)
  images=$(ls *.png)
  for image in $images; do
    date=$(echo $image | grep -Eo '_[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
    if ! [[ -z "${date// }" ]]; then
      mv "${image}" "${image/${date}/}"
    fi
  done

  # Remove date for quantile table (if quantile plot was produced)
  if ls PRSice_QUANTILES*.txt 1> /dev/null 2>&1; then
    table=$(ls PRSice_QUANTILES*.txt)
    date=$(echo $table | grep -Eo '_[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')
    if ! [[ -z "${date// }" ]]; then
      mv "${table}" "${table/${date}/}"
    fi
  fi
  '''
}



/*--------------------------
  Additional visualizations
----------------------------*/
 
process additional_plots {
  publishDir "${params.outdir}", mode: "copy"

  input:
  tuple file(pheno), file(cov) from transformed_target_pheno_for_plots_ch
  file prs from best_PRS_ch
  file metadata from pheno_metadata_ch

  output:
  file("*.png") into more_plots_ch

  script:
  """
  plot_prs_vs_cov.R --input_cov ${cov} --input_prs ${prs} --input_metadata ${metadata}
  plot_prs_density.R --input_pheno ${pheno} --input_prs ${prs}
  """

}



/*--------------------------
  Produce R Markdown report                          
----------------------------*/

process produce_report {
  publishDir params.outdir, mode: "copy"

  input:
  file plots from more_plots_ch
  file("*") from results_for_report_ch

  output:
  file ("MultiQC/multiqc_report.html") into reports

  script:
  if (params.quantile) {
    quantile_plot = "PRSice_QUANTILES_PLOT.png"
    quantile_table = "PRSice_QUANTILES.txt"
  } else {
    quantile_plot = "FALSE"
    quantile_table = "FALSE"
  }
  """
  cp /opt/bin/* .

  R -e "rmarkdown::render('prs_report.Rmd', params = list(barplot='PRSice_BARPLOT.png', highres.plot='PRSice_HIGH-RES_PLOT.png', density.plot='prs-density.png', quantile.plot='${quantile_plot}', quantile.table='${quantile_table}', prs.prsice='PRSice.prsice', prs.summary='PRSice.summary'))"
  mkdir MultiQC && mv prs_report.html MultiQC/multiqc_report.html

  """
}


