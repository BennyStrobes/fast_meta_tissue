import numpy as np 
import os
import sys
import pdb
import fast_meta_tissue as fmt

# Determine number of tests we are going to run
def extract_number_of_tests(trans_eqtl_file):
	f = open(trans_eqtl_file)
	head_count = 0  # Used to skip header
	count = 0  # Keep track of number of tests
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0: # skip header
			head_count = head_count + 1
			continue
		count = count + 1
	return count

# Determine test start number and test end number of this job
# For parrallelization purposes
def compute_test_start_and_end_number(job_number, num_jobs, num_tests):
	num_tests_per_job = num_tests/num_jobs + 1
	start_number = job_number*num_tests_per_job
	end_number =(job_number + 1)*num_tests_per_job - 1
	return start_number, end_number

trans_eqtl_file = sys.argv[1]
preprocess_metatissue_dir = sys.argv[2]
metatissue_output_dir = sys.argv[3]
java_compiler = sys.argv[4]
han_eskin_pvalue_table_file = sys.argv[5]
mt_pvalue_table_file = sys.argv[6]
metasoft_jar = sys.argv[7]
gemma_directory = sys.argv[8]
job_number = int(sys.argv[9])
num_jobs = int(sys.argv[10])

# Determine number of tests we are going to run
num_tests = extract_number_of_tests(trans_eqtl_file)

# Determine test start number and test end number of this job
# For parrallelization purposes
start_number, end_number = compute_test_start_and_end_number(job_number, num_jobs, num_tests)

heuristic = False

head_count = 0
counter = -1
f = open(trans_eqtl_file)
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	counter = counter + 1
	if counter >= start_number and counter <= end_number:
		snp_id = data[0]
		ensamble_id = data[1]
		gene_list_file = preprocess_metatissue_dir + 'gene_list_' + snp_id + '_' + ensamble_id + '.txt'
		genotype_file = preprocess_metatissue_dir + 'genotype_' + snp_id + '_' + ensamble_id + '.txt'
		tissue_info_file = preprocess_metatissue_dir + 'tissue_info_' + snp_id + '_' + ensamble_id + '.txt'
		print(snp_id + '\t' + ensamble_id)
		boolean = fmt.run_fast_meta_tissue(metatissue_output_dir, gene_list_file, genotype_file, tissue_info_file, heuristic, snp_id, ensamble_id, gemma_directory, java_compiler, metasoft_jar, han_eskin_pvalue_table_file, mt_pvalue_table_file)

