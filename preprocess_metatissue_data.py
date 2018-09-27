import numpy as np 
import os
import sys
import pdb
import gzip

# Extract list of tissues
def extract_gtex_tissues(file_name):
	arr = []
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		arr.append(line)
	return np.asarray(arr)


# First create object that maps from variants (found in eqtl_file) to:
## vector of covariate corrected dosage data
def initialize_variant_object(eqtl_list_file):
	variants = {}
	f = open(eqtl_list_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Snp id is first column
		snp_id = data[0]
		# Add variant to dictionary
		variants[snp_id] = []
	return variants

# create object that maps from gene names (found in eqtl file) to:
## Vector for each tissue of expression in that tissue
def initialize_expression_object(eqtl_list_file, tissues):
	genes = {}
	f = open(eqtl_list_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# gene id is second column
		gene_id = data[1]
		# Add gene to dictionary
		genes[gene_id] = {}
		# Add key to dictionary that will keep track of which tissues are observed for which gene
		genes[gene_id]['observed_tissues'] = []
		# For each tissue, add key that will keep track of expression measuremnets for that tissue
		for tissue in tissues:
			genes[gene_id][tissue] = []
	return genes

# get ordered list of donors from the genotype data files

def get_genotype_donors(file_name):
	f = gzip.open(file_name)
	head_count = 0
	for line in f:
		# All information need here can be found in first line of file
		if head_count == 0:
			head_count = head_count + 1
			line = line.rstrip()
			data = line.split()
			# Names of samples exclude first three elements
			donors = data[3:]
			continue
		break
	f.close()
	# Make sure donors are in alphabetical order!!
	if np.array_equal(donors,sorted(donors)) == False:
		print('fatal assumption made in processing variant data. Need to quit')
		pdb.set_trace()
	return np.asarray(donors)

# Get ordered lists of donors (one for each tissue) from the expression data files
def get_expression_donors(expression_dir, tissues):
	donors = {}
	# Loop through tisues
	for tissue in tissues:
		# Expression file for this tissue
		file_name = expression_dir + tissue + '.v8.normalized_expression_adjusted_by_covariates.txt.gz'
		head_count = 0
		f = gzip.open(file_name)
		for line in f:
			# All information we need can be found in the first line of this file
			if head_count == 0:
				head_count = head_count + 1
				line = line.rstrip()
				data = line.split()
				tissue_specific_donors = data[4:]
				continue
			break
		f.close()
		# Make sure donor names are in alphabetical order
		if np.array_equal(tissue_specific_donors,sorted(tissue_specific_donors)) == False:
			print('fatal assumption made in processing expression data. Need to quit')
			pdb.set_trace()
		# Add tissue specific donors to dictionary
		donors[tissue] = np.asarray(tissue_specific_donors)
	return donors

# Fill in variant data for current chromosome
def fill_in_variant_data_one_chromosome(variants, genotype_file):
	f = gzip.open(genotype_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_id = data[2]
		# Only consider variants that are in variants dictionary. Otherwise, ignore
		if snp_id not in variants:
			continue
		genotype_vec = data[3:]
		variants[snp_id] = genotype_vec
	return variants

# Loop through all variants of interest and make sure we collected data (error checking)
def check_to_make_sure_all_variants_were_filled_in(variants):
	for variant in variants.keys():
		if len(variants[variant]) == 0:
			print('variant not filled in! Must hault.')
			pdb.set_trace()
	return


# Fill in variant data
def fill_in_variant_data(variants, genotype_dir):
	# Do seperately for each chromosome
	for chrom_num in range(1,23):
		genotype_file = genotype_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + str(chrom_num) + '_adjusted_dosage.txt.gz'
		# Fill in variant data for current chromosome
		variants = fill_in_variant_data_one_chromosome(variants, genotype_file)
	# Loop through all variants of interest and make sure we collected data (error checking)
	check_to_make_sure_all_variants_were_filled_in(variants)
	return variants

# Fill in expression data for this one tissue
def fill_in_expression_data_for_one_tissue(expression_file, genes, tissue):
	f = gzip.open(expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[3]
		# Ignore genes that are not in our list of $genes
		if ensamble_id not in genes:
			continue
		# Extract vector of gene expression measurements
		expression_vec = data[4:]
		# Add expression vec to genes dictionary
		genes[ensamble_id][tissue] = expression_vec
		# Add tissue to list of observed tissues for this gene
		genes[ensamble_id]['observed_tissues'].append(tissue)
	return genes

# Fill in expression data
def fill_in_expression_data(genes, expression_dir, tissues):
	# Loop through tissues
	for tissue in tissues:
		# Expression file for this tissues
		expression_file = expression_dir + tissue + '.v8.normalized_expression_adjusted_by_covariates.txt.gz'
		# Fill in expression data for this one tissue
		genes = fill_in_expression_data_for_one_tissue(expression_file, genes, tissue)
	return genes

# Extract list of variant gene pairs
def extract_list_of_variant_gene_pairs(eqtl_list_file):
	f = open(eqtl_list_file)
	variant_gene_pairs = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# extract variant and ensamble ids
		snp_id = data[0]
		ensamble_id = data[1]
		# Create variant-gene pair test name with hyphen merge
		test_name = snp_id + '-' + ensamble_id
		# Add test_name to array
		variant_gene_pairs.append(test_name)
	return np.asarray(variant_gene_pairs)

# Print genotype data to output file
def print_genotype_data(genotype_vector, genotype_output_file):
	t = open(genotype_output_file, 'w')
	t.write('\t'.join(genotype_vector) + '\n')
	t.close()
	return

# Print expression data for each tissue
def print_expression_data(tissue_expression_output_file, tissue_specific_donors, tissue_specific_expression):
	# Check to make sure lenth's are equivalent
	if len(tissue_specific_expression) != len(tissue_specific_donors):
		print('Length mismatch! Must stop!!')
	# Open output file handle
	t = open(tissue_expression_output_file, 'w')
	# Write header to output file handle
	t.write('\t'.join(tissue_specific_donors) + '\n')
	# Write expression vec to output file handle
	t.write('\t'.join(tissue_specific_expression) + '\n')
	t.close()

def make_global_tissue_info_mat(genotype_donors, expression_donors, tissues):
	tissue_info_mat = np.zeros((len(genotype_donors), len(tissues)))
	print(len(genotype_donors))
	for i, genotype_donor in enumerate(genotype_donors):
		for j, tissue in enumerate(tissues):
			if genotype_donor in expression_donors[tissue]:
				tissue_info_mat[i,j] = 1.0
	return tissue_info_mat

# Make tissue info file for this variant gene pair
def print_tissue_info_file(tissue_info_file, tissue_info_mat, all_tissues, observed_tissues, genotype_donors):
	# Extract indices of all_tissues vector whose elements are in observed_tissues
	valid_column_indices = []
	for i, tissue in enumerate(all_tissues):
		if tissue in observed_tissues:
			valid_column_indices.append(i)
	# Subset tissue_info_mat to only include columns (tissues) that are observed
	filtered_tissue_info_mat = tissue_info_mat[:, valid_column_indices]
	# Convert tissue info matrix to be of type string
	filtered_tissue_info_mat = filtered_tissue_info_mat.astype(int).astype(str)
	# Add row labels
	filtered_tissue_info_mat = np.hstack((np.transpose(np.asmatrix(genotype_donors)), filtered_tissue_info_mat))
	# Add column labels
	header = np.asmatrix(['#TISSUE'] + observed_tissues)
	filtered_tissue_info_mat = np.vstack((header,filtered_tissue_info_mat))
	# Save matrix to file
	np.savetxt(tissue_info_file, filtered_tissue_info_mat, fmt="%s",delimiter='\t')
	return

def print_helper(variant_gene_pairs, tissues, variants, genes, genotype_donors, expression_donors, preprocess_metatissue_dir):
	tissue_info_mat = make_global_tissue_info_mat(genotype_donors, expression_donors, tissues)
	# Loop through variant gene pairs
	for variant_gene_pair in variant_gene_pairs:
		print(variant_gene_pair)
		# Variant Id
		snp_id = variant_gene_pair.split('-')[0]
		# Ensamble id
		ensamble_id = variant_gene_pair.split('-')[1]
		# Extract list of observed tissues for this gene
		observed_tissues = genes[ensamble_id]['observed_tissues']
		# Print genotype data to output file
		genotype_output_file = preprocess_metatissue_dir + 'genotype_' + snp_id + '_' + ensamble_id + '.txt'
		print_genotype_data(variants[snp_id], genotype_output_file)

		# Open filehandle for gene_list (each line of file is file_name of expression for ONE tissue)
		gene_list_file_handle = open(preprocess_metatissue_dir + 'gene_list_' + snp_id + '_' + ensamble_id + '.txt', 'w')
		# Print expression data for each tissue
		for tissue in observed_tissues:
			# Print expression for this tissue
			tissue_expression_output_file = preprocess_metatissue_dir + 'expression_' + tissue + '_' + snp_id + '_' + ensamble_id + '.txt'
			print_expression_data(tissue_expression_output_file, expression_donors[tissue], genes[ensamble_id][tissue])
			# Add tissue specific expression file to gene_list
			gene_list_file_handle.write(tissue_expression_output_file + '\n')
		gene_list_file_handle.close()
	
		# Make tissue info file for this variant gene pair
		tissue_info_file = preprocess_metatissue_dir + 'tissue_info_' + snp_id + '_' + ensamble_id + '.txt'
		print_tissue_info_file(tissue_info_file, tissue_info_mat, tissues, observed_tissues, genotype_donors)





eqtl_list_file = sys.argv[1]
tissue_list_file = sys.argv[2]
genotype_dir = sys.argv[3]
expression_dir = sys.argv[4]
preprocess_metatissue_dir = sys.argv[5]


# Extract list of tissues
tissues = extract_gtex_tissues(tissue_list_file)

# Extract list of variant gene pairs
variant_gene_pairs = extract_list_of_variant_gene_pairs(eqtl_list_file)

# First create object that maps from variants (found in eqtl_file) to:
## vector of covariate corrected dosage data
variants = initialize_variant_object(eqtl_list_file)

# create object that maps from gene names (found in eqtl file) to:
## Vector for each tissue of expression in that tissue
genes = initialize_expression_object(eqtl_list_file, tissues)


# get ordered list of donors from the genotype data files
genotype_donors = get_genotype_donors(genotype_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr1_adjusted_dosage.txt.gz')

# Get ordered lists of donors (one for each tissue) from the expression data files
expression_donors = get_expression_donors(expression_dir, tissues)

# Fill in variant data
variants = fill_in_variant_data(variants, genotype_dir)

# Fill in expression data
genes = fill_in_expression_data(genes, expression_dir, tissues)


print_helper(variant_gene_pairs, tissues, variants, genes, genotype_donors, expression_donors, preprocess_metatissue_dir)