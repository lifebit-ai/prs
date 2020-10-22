/*------------------------------------------------------------------------
  Detecting target and base builds and updating the base build if need be 
--------------------------------------------------------------------------*/

header_ch = Channel
      .fromPath(params.header, checkIfExists: true)
      .ifEmpty { exit 1, "Header file required for detecting builds not found: ${params.header}" }



process detect_and_update_build {
    publishDir "${params.outdir}/transformed_PRSice_inputs", mode: "copy"
    
    input:
    tuple val(name), file("*") from target_plink_build_ch
    file base from transformed_base_ch
    file header from header_ch
    
    output:
    file("matching.base.data") into matching_base_build_ch
    
    shell:
    '''
    # Step 1 - Combine all input *bim files into 1 file and transform to a 23andme format so that Python package "snps" can read it
    
    cat *bim \
    | cut -f 1,2,4,5,6 \
    | awk '{print $2"\t"$1"\t"$3"\t"$4 $5}' > all_bims.tsv
    cat !{header} all_bims.tsv > all_bims-transformed.tsv
    
    # Step 2 - Transform the base.data file to a 23andme format so that Python package "snps" can read it

    cat base.data \
    | sed '1d' \
    | cut -f 1,2,3,4,5 -d " " \
    | awk '{print $3"\t"$1"\t"$2"\t"$4 $5}' > base.tsv
    cat !{header} base.tsv > base-transformed.tsv

    # Step 3 - Use Python package "snsp" to detect build of target and base and update base build if need be
    
    detect_and_update_build.py --input_target all_bims-transformed.tsv --input_base base-transformed.tsv

    # Step 4 - If step 3 produced a file, use it to update the build of the base

    if ls new_base_coordinates.txt 1> /dev/null 2>&1
    then
      update_base_build.R --input_base base.data --input_coordinates new_base_coordinates.txt
    else
      mv base.data matching.base.data
    fi
    '''
}
