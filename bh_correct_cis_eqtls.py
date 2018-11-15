import numpy as np 
import os
import sys
import pdb
import random



# Extract array of tissue names from file
def get_tissues(tissue_list_file):
	tissues = []
	f = open(tissue_list_file)
	for line in f:
		line = line.rstrip()
		tissues.append(line)
	return np.asarray(tissues)

def bh_correction(cis_eqtl_file, bh_corrected_cis_eqtl_file, fdr_thresh, used_variants, used_genes):
	# Initialize array to keep track of cis-eqtls and their pvalues
	cis_eqtls = []
	# Used to skip header
	head_count = 0
	# Used to keep track of total number of tests
	num_tests = 0

	# Stream cis-eqtl file
	f = open(cis_eqtl_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent info for this variant-gene pair
		ensamble_id = data[0]
		snp_id = data[1]
		# Skip variant-gene pairs that were not used in the trans analysis
		if ensamble_id not in used_genes or snp_id not in used_variants:
			continue
		tss_distance = data[2]
		maf = data[5]
		pvalue = float(data[6])
		# Keep track of number of tests
		num_tests = num_tests + 1
		# Can throw these out won't be sig
		if pvalue > float(fdr_thresh):
			continue
		# Add relevent info to array
		cis_eqtls.append((pvalue, ensamble_id, snp_id, tss_distance, maf))
	sorted_cis_eqtls = sorted(cis_eqtls, key=lambda tup: tup[0])

	iteration = 0
	converged = False
	t = open(bh_corrected_cis_eqtl_file, 'w')
	t.write('snp_id\tensamble_id\tss_distance\tmaf\tpvalue\tfdr\n')
	for variant_gene_pair_tuple in sorted_cis_eqtls:
		iteration = iteration + 1
		pvalue = variant_gene_pair_tuple[0]
		ensamble_id = variant_gene_pair_tuple[1]
		snp_id = variant_gene_pair_tuple[2]
		tss_distance = variant_gene_pair_tuple[3]
		maf = variant_gene_pair_tuple[4]
		fdr = (pvalue*num_tests)/iteration
		if fdr > fdr_thresh:
			converged = True
			break
		t.write(snp_id + '\t' + ensamble_id + '\t' + tss_distance + '\t' + maf + '\t' + str(pvalue) + '\t' + str(fdr) + '\n')
	t.close()
	if converged == False:
		print('EROROROR')
		print(cis_eqtl_file)
		pdb.set_trace()
	return

# Open output file handles and write header
def open_handle(file_name):
	t = open(file_name,'w')
	t.write('snp_id\tensamble_id\ttss_distance\tmaf\tpvalue\tfdr\n')
	return t


def get_dictionary_list_of_test_names(bh_corrected_cis_eqtl_file):
	# Mapping from gene name to array of all variants in that gene
	gene_to_variants = {}
	# Mapping from gene to top cis variant
	gene_to_top_cis_variant = {}

	# Open eqtl file
	f = open(bh_corrected_cis_eqtl_file)
	head_count = 0  # to skip header
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		head_count = head_count + 1
		# Extract relevent fields from line
		snp_id = data[0]
		ensamble_id = data[1]
		pvalue = float(data[-2])		

		# Add variant to list of variants for this gene
		if ensamble_id not in gene_to_variants:
			gene_to_variants[ensamble_id] = []
		gene_to_variants[ensamble_id].append(snp_id)
		# Update top variant per gene
		if ensamble_id not in gene_to_top_cis_variant:
			gene_to_top_cis_variant[ensamble_id] = (snp_id, pvalue)
		else:
			old_pvalue = gene_to_top_cis_variant[ensamble_id][1]
			if pvalue < old_pvalue:
				gene_to_top_cis_variant[ensamble_id] = (snp_id, pvalue)
	f.close()
	egenes_and_random = {}
	egenes_and_top = {}
	# get list of egenes and random
	for egene in gene_to_variants.keys():
		variants = gene_to_variants[egene]
		randomly_selected_variant = random.choice(variants)
		egenes_and_random[randomly_selected_variant + '_' + egene] = 1
	# Get list of egenes and top variant
	for egene in gene_to_top_cis_variant.keys():
		variant = gene_to_top_cis_variant[egene][0]
		egenes_and_top[variant + '_' + egene] = 1

	return egenes_and_random, egenes_and_top


def print_egene_files(t_cis, t_top_cis, bh_corrected_cis_eqtl_file):
	egenes_and_random, egenes_and_top = get_dictionary_list_of_test_names(bh_corrected_cis_eqtl_file)
	f = open(bh_corrected_cis_eqtl_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		test_name = data[0] + '_' + data[1]
		if test_name in egenes_and_random:
			t_cis.write(line + '\n')
		if test_name in egenes_and_top:
			t_top_cis.write(line + '\n')
	f.close()
	t_cis.close()
	t_top_cis.close()
	return

def extract_list_of_variants_to_be_used_in_this_analysis(trans_variant_list_file):
	f = open(trans_variant_list_file)
	variants = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variants[data[0]] = 1
	return variants

def extract_list_of_genes_to_be_used_in_this_analysis(trans_gene_list_file):
	f = open(trans_gene_list_file)
	genes = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes[data[0]] = 1
	return genes

tissue = sys.argv[1]
raw_cis_eqtl_results_dir = sys.argv[2]
cis_eqtl_bh_corrected_dir = sys.argv[3]
fdr_thresh = float(sys.argv[4])
trans_variant_list_file = sys.argv[5]
trans_gene_list_file = sys.argv[6]

print(tissue)
print(fdr_thresh)

used_variants = extract_list_of_variants_to_be_used_in_this_analysis(trans_variant_list_file)
used_genes = extract_list_of_genes_to_be_used_in_this_analysis(trans_gene_list_file)

# Input file for tissue
cis_eqtl_file = raw_cis_eqtl_results_dir + tissue + '.allpairs.txt'
# Output file for tissue
bh_corrected_cis_eqtl_file = cis_eqtl_bh_corrected_dir + tissue + '_bh_corrected_' + str(fdr_thresh) + '.txt'
	
# Run bh correction on this tissues cis-eqtls using variant-gene pair level fdr correction (not gene level)
bh_correction(cis_eqtl_file, bh_corrected_cis_eqtl_file, fdr_thresh, used_variants, used_genes)



# Extract two list of eqtls from these files:
## 1 for all eGenes and a randomly selected variant of theirs that passes FDR <= FDR thresh
## 2. For all egenes and their top associated variant
# Open output file handles and write header
t_cis = open_handle(cis_eqtl_bh_corrected_dir + tissue + '_egenes_and_random_variant_fdr_' + str(fdr_thresh) + '.txt')
t_top_cis = open_handle(cis_eqtl_bh_corrected_dir + tissue + '_egenes_and_top_variant_fdr_' + str(fdr_thresh) + '.txt')
print_egene_files(t_cis, t_top_cis, bh_corrected_cis_eqtl_file)
