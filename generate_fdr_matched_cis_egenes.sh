#!/bin/bash -l

#SBATCH
#SBATCH --time=4:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1

trans_eqtl_file="$1"
tissue_list_file="$2"
cis_eqtl_bh_corrected_dir="$3"
cis_eqtl_file="$4"
top_cis_eqtl_file="$5"




python generate_fdr_matched_cis_egenes.py $trans_eqtl_file $tissue_list_file $cis_eqtl_bh_corrected_dir $cis_eqtl_file $top_cis_eqtl_file