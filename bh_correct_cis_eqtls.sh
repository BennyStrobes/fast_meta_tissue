#!/bin/bash -l

#SBATCH
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=40GB

tissue_name="$1"
raw_cis_eqtl_results_dir="$2"
cis_eqtl_bh_corrected_dir="$3"
fdr_thresh="$4"
trans_variant_list_file="$5"
trans_gene_list_file="$6"
trans_eqtl_file="$7"



if false; then
python bh_correct_cis_eqtls.py $tissue_name $raw_cis_eqtl_results_dir $cis_eqtl_bh_corrected_dir $fdr_thresh $trans_variant_list_file $trans_gene_list_file
fi



python find_cis_trans_overlaps.py $tissue_name $cis_eqtl_bh_corrected_dir $fdr_thresh $trans_eqtl_file