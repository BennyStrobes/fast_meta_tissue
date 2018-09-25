
trans_eqtl_file="$1"
preprocess_metatissue_dir="$2"
metatissue_output_dir="$3"
java_compiler="$4"
han_eskin_pvalue_table_file="$5"
mt_pvalue_table_file="$6"
metasoft_jar="$7"
gemma_directory="$8"

python run_gemma_tissue.py $trans_eqtl_file $preprocess_metatissue_dir $metatissue_output_dir $java_compiler $han_eskin_pvalue_table_file $mt_pvalue_table_file $metasoft_jar $gemma_directory