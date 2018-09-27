#!/bin/bash -l

#SBATCH
#SBATCH --time=2:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


eqtl_file="$1"
organized_metatissue_output_dir="$2"
metatissue_output_dir="$3"
preprocess_metatissue_dir="$4"
tissue_list_file="$5"


python organize_metatissue_results.py $eqtl_file $organized_metatissue_output_dir $metatissue_output_dir $preprocess_metatissue_dir $tissue_list_file