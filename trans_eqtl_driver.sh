#!/bin/bash -l

#SBATCH
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB


################################################
## INPUT DATA
#################################################
# File containing list of variant-gene pairs we wish to test
# Each line is a variant gene pair
# Column 1 is variant id
# Column 2 is ensamble id
trans_eqtl_file="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/cis_practice_input.txt"

# File containing all gtex v8 tissue names (each line is a tissue)
tissue_list_file="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/tissue_names.txt"

# Directory containing all genotype data for v8 (where covariates have already been regressed out)
# One file per chromosome
# Files of format GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr${chrom_num}_adjusted_dosage.txt.gz
genotype_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/"

# Directory containing expression file for each of gtex tissues
# Files of format ${tissue_name}.v8.normalized_expression_adjusted_by_covariates.txt.gz
expression_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/gtex_v8_gex_adjusted/"

# Location of java compiler
java_compiler="/software/cmshared/apps/java/jdk1.8.0_112/bin/java"

# Set file from metatissue package
han_eskin_pvalue_table_file="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/HanEskinPvalueTable.txt"

# Set file from metatissue package
mt_pvalue_table_file="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/MTPvalueTable.txt.gz"

# Metasoft jar file
metasoft_jar="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/Metasoft.jar"

# Directory containing gemma executable
gemma_directory="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/"









################################################
## Output DATA
#################################################

# Directory containing files used as input to metatissue
# Created in preprocess_metatissue_data.sh
preprocess_metatissue_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/metatissue_input/"


# Directory containing metatissue output files
metatissue_output_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/metatissue_output/"

# Directory containing organized (merged) metatissue results
organized_metatissue_output_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/organized_metatissue_results/"






################################################
## Run analysis
#################################################
sh preprocess_metatissue_data.sh $trans_eqtl_file $tissue_list_file $genotype_dir $expression_dir $preprocess_metatissue_dir

sh run_gemma_tissue.sh $trans_eqtl_file $preprocess_metatissue_dir $metatissue_output_dir $java_compiler $han_eskin_pvalue_table_file $mt_pvalue_table_file $metasoft_jar $gemma_directory


sh organize_metatissue_results.sh $trans_eqtl_file $organized_metatissue_output_dir $metatissue_output_dir $preprocess_metatissue_dir $tissue_list_file