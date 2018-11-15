#!/bin/bash
#SBATCH --time=24:00:00 --partition=broadwl --mem=20GB



################################################
## INPUT DATA
#################################################
# File containing list of variant-gene pairs we wish to test
# Each line is a variant gene pair
# Column 1 is variant id
# Column 2 is ensamble id
#cis_practice_input.txt
trans_eqtl_file="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/protein_coding_trans_hits_fdr_.5.txt"
trans_egene_file="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/input_data/protein_coding_trans_egenes_fdr_.5.txt"

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

# Directory containing all raw cis-eqtl results
# File for each tissue of format $tissue_name".allpairs.txt"
raw_cis_eqtl_results_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/"


# Directory containing list of variants used in trans analysis (seperate file per tissue)
trans_variants_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/input_data/trans_variant_list/"

# Directory containing list of genes used in trans analysis (seperate file per tissue)
trans_genes_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/input_data/trans_gene_list/"






################################################
## Output DATA
#################################################

# Directory containing files used as input to metatissue for trans eqtls
# Created in preprocess_metatissue_data.sh
preprocess_trans_metatissue_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/trans_metatissue_input/"

# Directory containing metatissue output files for tans eqtls
trans_metatissue_output_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/trans_metatissue_output/"

# Directory containing organized (merged) metatissue results
organized_metatissue_output_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/organized_metatissue_results/"

# Directory containing cis-eQTLs based on Benjamini-hochberg correction on ALL Variant gene pairs
cis_eqtl_bh_corrected_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/cis_eqtl_bh_corrected/"

# Directory containing files used as input to metatissue for matched cis-eqtls
preprocess_cis_metatissue_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/cis_metatissue_input/"

# Directory containing metatissue output files for tans eqtls
cis_metatissue_output_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/cis_metatissue_output/"

# Directory containing metatissue output files for tans eqtls
top_cis_metatissue_output_dir="/work-zfs/abattle4/bstrober/gtex_v8_eqtls/tissue_specificity/metatissue/top_cis_metatissue_output/"




################################################
## Create list of cis-eqtls after BH correction across all variant gene pairs (not at gene level) to make comparible with trans results
# Do for two FDR thresholds
#################################################
fdr_thresh=".25"
if false; then
while read tissue_name; do
	trans_gene_list_file=$trans_genes_dir$tissue_name"_gene_anno.tsv"
	trans_variant_list_file=$trans_variants_dir$tissue_name"_variant_list.tsv"
	sbatch bh_correct_cis_eqtls.sh $tissue_name $raw_cis_eqtl_results_dir $cis_eqtl_bh_corrected_dir $fdr_thresh $trans_variant_list_file $trans_gene_list_file $trans_egene_file
done <$tissue_list_file
fi

tissue_name="Adipose_Subcutaneous"
trans_gene_list_file=$trans_genes_dir$tissue_name"_gene_anno.tsv"
trans_variant_list_file=$trans_variants_dir$tissue_name"_variant_list.tsv"
sh bh_correct_cis_eqtls.sh $tissue_name $raw_cis_eqtl_results_dir $cis_eqtl_bh_corrected_dir $fdr_thresh $trans_variant_list_file $trans_gene_list_file $trans_egene_file

fdr_thresh=".5"
if false; then
while read tissue_name; do
	trans_gene_list_file=$trans_genes_dir$tissue_name"_gene_anno.tsv"
	trans_variant_list_file=$trans_variants_dir$tissue_name"_variant_list.tsv"
	sbatch bh_correct_cis_eqtls.sh $tissue_name $raw_cis_eqtl_results_dir $cis_eqtl_bh_corrected_dir $fdr_thresh $trans_variant_list_file $trans_gene_list_file $trans_egene_file
done <$tissue_list_file
fi











################################################
## Run metatissue (gemma implementation) analysis for Trans eQTLS
#################################################

#####################################
# Part 1: Generate metatissue input files for trans-eqtls
if false; then
sh preprocess_metatissue_data.sh $trans_eqtl_file $tissue_list_file $genotype_dir $expression_dir $preprocess_trans_metatissue_dir
fi

#####################################
# Part 2: Run fast metatissue 
num_jobs="400"
#for job_number in $(seq 0 $(($num_jobs-1))); do 
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	echo $job_number
	sbatch run_gemma_tissue.sh $trans_eqtl_file $preprocess_trans_metatissue_dir $trans_metatissue_output_dir $java_compiler $han_eskin_pvalue_table_file $mt_pvalue_table_file $metasoft_jar $gemma_directory $job_number $num_jobs
done
fi

#####################################
# Part 3: Merge results from one metatissue run into a single output file for that run

organized_trans_eqtl_metatissue_results_file=$organized_metatissue_output_dir"gtex_v8_trans_eqtl_gemmatissue_fdr_.5.txt"
if false; then
sh organize_metatissue_results.sh $trans_eqtl_file $organized_trans_eqtl_metatissue_results_file $trans_metatissue_output_dir $preprocess_trans_metatissue_dir $tissue_list_file 
fi

















################################################
## Run metatissue (gemma implementation) analysis for FDR-matched Cis eQTLS
#################################################


cis_eqtl_file=$preprocess_cis_metatissue_dir"protein_coding_cis_egenes_fdr_.5.txt"
top_cis_eqtl_file=$preprocess_cis_metatissue_dir"protein_coding_top_cis_egenes_fdr_.5.txt"


#####################################
# Part 1: Generate two randomly selected list of cis-eqtl genes:
## 1 for just random variant-eGene pairs
## 2. For random eGenes (and their top cis associated variant)
if false; then
sh generate_fdr_matched_cis_egenes.sh $trans_egene_file $tissue_list_file $cis_eqtl_bh_corrected_dir $cis_eqtl_file $top_cis_eqtl_file
fi

#####################################
# Part 2: Generate metatissue input files for cis-eqtls

# for randomly selected variants
if false; then
sbatch preprocess_metatissue_data.sh $cis_eqtl_file $tissue_list_file $genotype_dir $expression_dir $preprocess_cis_metatissue_dir
# For top variants
sbatch preprocess_metatissue_data.sh $top_cis_eqtl_file $tissue_list_file $genotype_dir $expression_dir $preprocess_cis_metatissue_dir
fi

#####################################
# Part 3: Run fast metatissue 
# For randomly matched cis-eqtls
num_jobs="100"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	echo $job_number
	sbatch run_gemma_tissue.sh $cis_eqtl_file $preprocess_cis_metatissue_dir $cis_metatissue_output_dir $java_compiler $han_eskin_pvalue_table_file $mt_pvalue_table_file $metasoft_jar $gemma_directory $job_number $num_jobs
done
fi
# And randomly matched top cis-eqtls
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	echo $job_number
	sbatch run_gemma_tissue.sh $top_cis_eqtl_file $preprocess_cis_metatissue_dir $top_cis_metatissue_output_dir $java_compiler $han_eskin_pvalue_table_file $mt_pvalue_table_file $metasoft_jar $gemma_directory $job_number $num_jobs
done
fi



#####################################
# Part 4: Merge results from one metatissue run into a single output file for that run
if false; then
organized_cis_eqtl_metatissue_results_file=$organized_metatissue_output_dir"gtex_v8_cis_egenes_random_variant_eqtl_gemmatissue_fdr_.5.txt"
sbatch organize_metatissue_results.sh $cis_eqtl_file $organized_cis_eqtl_metatissue_results_file $cis_metatissue_output_dir $preprocess_cis_metatissue_dir $tissue_list_file 


organized_cis_eqtl_metatissue_results_file=$organized_metatissue_output_dir"gtex_v8_cis_egenes_top_variant_eqtl_gemmatissue_fdr_.5.txt"
sbatch organize_metatissue_results.sh $top_cis_eqtl_file $organized_cis_eqtl_metatissue_results_file $top_cis_metatissue_output_dir $preprocess_cis_metatissue_dir $tissue_list_file 
fi




















################################################
## Plotting (temporarily running on RCC)
#################################################

# Metatissue results using cis-eqtl results from YoSon
all_cis_eqtl_file="/project2/gilad/bstrober/metatissue/input_data/gemmatissue_output_w_gtex_v8_300k.txt"
# Metatissue results using trans-eqtl results
all_trans_eqtl_file="/project2/gilad/bstrober/metatissue/input_data/gtex_v8_trans_eqtl_gemmatissue_fdr_.5.txt"
# Metatissue results using trans-egenes and their top associated variant
trans_egene_file="/project2/gilad/bstrober/metatissue/input_data/gtex_v8_trans_egenes_gemmatissue_fdr_.5.txt"
# Metatissue results using random cis-egenes (FDR <= .5) a random accociated variant
fdr_matched_cis_egene_file="/project2/gilad/bstrober/metatissue/input_data/gtex_v8_cis_egenes_random_variant_eqtl_gemmatissue_fdr_.5.txt"
# Metatissue results using random cis-egenes (FDR <= .5) and their top associated variant
fdr_matched_top_cis_egene_file="/project2/gilad/bstrober/metatissue/input_data/gtex_v8_cis_egenes_top_variant_eqtl_gemmatissue_fdr_.5.txt"
# tissue names file
tissue_names_file="/project2/gilad/bstrober/metatissue/input_data/tissue_names.txt"
# File containing gtex tissue colors
gtex_tissue_colors_file="/project2/gilad/bstrober/metatissue/input_data/gtex_colors.txt"

output_dir="/project2/gilad/bstrober/metatissue/output_data/"

if false; then
Rscript make_tissue_sharing_heatmap.R $all_cis_eqtl_file $all_trans_eqtl_file $tissue_names_file $gtex_tissue_colors_file $output_dir
fi
if false; then
Rscript make_mvalue_distribution_histogram.R $trans_egene_file $fdr_matched_cis_egene_file $fdr_matched_top_cis_egene_file $output_dir
fi



