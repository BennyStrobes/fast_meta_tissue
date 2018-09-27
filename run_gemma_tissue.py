import numpy as np 
import os
import sys
import pdb
import fast_meta_tissue as fmt





trans_eqtl_file = sys.argv[1]
preprocess_metatissue_dir = sys.argv[2]
metatissue_output_dir = sys.argv[3]
java_compiler = sys.argv[4]
han_eskin_pvalue_table_file = sys.argv[5]
mt_pvalue_table_file = sys.argv[6]
metasoft_jar = sys.argv[7]
gemma_directory = sys.argv[8]

heuristic = False

head_count = 0
f = open(trans_eqtl_file)
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	snp_id = data[0]
	ensamble_id = data[1]
	gene_list_file = preprocess_metatissue_dir + 'gene_list_' + snp_id + '_' + ensamble_id + '.txt'
	genotype_file = preprocess_metatissue_dir + 'genotype_' + snp_id + '_' + ensamble_id + '.txt'
	tissue_info_file = preprocess_metatissue_dir + 'tissue_info_' + snp_id + '_' + ensamble_id + '.txt'

	boolean = fmt.run_fast_meta_tissue(metatissue_output_dir, gene_list_file, genotype_file, tissue_info_file, heuristic, snp_id, ensamble_id, gemma_directory, java_compiler, metasoft_jar, han_eskin_pvalue_table_file, mt_pvalue_table_file)

