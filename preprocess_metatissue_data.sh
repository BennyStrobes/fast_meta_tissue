
eqtl_list_file="$1"
tissue_list_file="$2"
genotype_dir="$3"
expression_dir="$4"
preprocess_metatissue_dir="$5"

python preprocess_metatissue_data.py $eqtl_list_file $tissue_list_file $genotype_dir $expression_dir $preprocess_metatissue_dir