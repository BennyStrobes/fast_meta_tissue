import numpy as np
import os
import sys
import pdb
import gzip

# Extract vector of tissue names from tissue_list_file input file
def extract_tissue_names(tissue_list_file):
	arr = []
	f = open(tissue_list_file)
	for line in f:
		line = line.rstrip()
		arr.append(line)
	return np.asarray(arr)


# Print header line only to output file
def print_header_to_output_file(t, tissue_names):
	t.write('RSID\tNUM_STUDY\tPVALUE_FE\tPVALUE_RE2\tI_SQUARE\tPVALUE_Q')
	for tissue_name in tissue_names:
		t.write('\t' + tissue_name + '_mvalue_stderr')
	for tissue_name in tissue_names:
		t.write('\t' + tissue_name + '_mvalue')
	for tissue_name in tissue_names:
		t.write('\t' + tissue_name + '_beta')
		t.write('\t' + tissue_name + '_stderr')
	t.write('\n')
	return t

# Extract dictionary list of all tissues observed for this variant gene pair from tissue info file
def extract_observed_tissues_from_tissue_info_file(tissue_info_file):
	f = open(tissue_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Observed tissue names are all in header.
		if head_count == 0:
			head_count = head_count + 1
			# Create array of observed tissue names
			tissue_arr = data[1:]
			# Convert from observed tissue array to dictionary
			tissue_dicti = {}
			for tissue in tissue_arr:
				tissue_dicti[tissue] = 1
			continue
		# Can break out of loop after header
		break
	return tissue_dicti



# Extract relevent fields from metasoft output file
# This includes:
#####1. PVALUE_FE
#####2. PVALUE_RE2
#####3. I_SQUARE
#####4. PVALUE_Q
#####5. vector of metatissue standard errors
#####6. vector of metatissue mvalues
def extract_relevent_fields_from_metasoft_output(metasoft_output_file, tissue_names, observed_tissues):
	num_observed_tissues = len(observed_tissues)
	# Open file handle
	f = gzip.open(metasoft_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:  # Skip header
			head_count = head_count + 1
			continue
		# extract relvent fields
		pvalue_fe = data[2]
		pvalue_re2 = data[3]
		i_square = data[4]
		pvalue_q = data[5]
		observed_mvalue_stderr_vector = data[6:(6+num_observed_tissues)]
		observed_mvalue_vector = data[(6+num_observed_tissues):]
	f.close()

	# Error checking
	if len(observed_mvalue_stderr_vector) != num_observed_tissues or len(observed_mvalue_vector) != num_observed_tissues:
		print('ASSUMPTION ERROR IN MERGING FILES. MUST STOP')
		pdb.set_trace()
	
	# Add NA's to mvalue vector and mvalue_stderr vector
	mvalue_vector = fill_in_vector_with_nas_for_missing_tissues(observed_mvalue_vector, tissue_names, observed_tissues)
	mvalue_stderr_vector = fill_in_vector_with_nas_for_missing_tissues(observed_mvalue_stderr_vector, tissue_names, observed_tissues)
	return pvalue_fe, pvalue_re2, i_square, pvalue_q, mvalue_stderr_vector, mvalue_vector


# Add NA's to vector if tissue is missing
def fill_in_vector_with_nas_for_missing_tissues(observed_mvalue_vector, tissue_names, observed_tissues):
	full_vector = []
	counter = 0
	# Loop through tissuses
	for tissue in tissue_names:
		# Tissue is observed for this gene
		if tissue in observed_tissues:
			full_vector.append(observed_mvalue_vector[counter])
			counter = counter + 1
		# Tissue is not observed for this gene
		else:
			full_vector.append('NA')
	return np.asarray(full_vector)


# Extract relevent fields from metatissue beta file
# This includes
#####1. Beta vector
#####2. stderr vector
def extract_relevent_fields_from_metatissue_beta_file(metatissue_beta_file, tissue_names, observed_tissues):
	num_observed_tissues = len(observed_tissues)

	# Loop through metatissue beta file
	f = gzip.open(metatissue_beta_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		observed_beta_vec = data[1:][::2]
		observed_stderr_vec = data[2:][::2]
	f.close()

	# Error checking
	if len(observed_beta_vec) != num_observed_tissues or len(observed_stderr_vec) != num_observed_tissues:
		print('ASSUMPTION ERROR IN MERGING FILES. MUST STOP')
		pdb.set_trace()

	# Add NA's to mvalue vector and mvalue_stderr vector
	beta_vec = fill_in_vector_with_nas_for_missing_tissues(observed_beta_vec, tissue_names, observed_tissues)
	stderr_vec = fill_in_vector_with_nas_for_missing_tissues(observed_stderr_vec, tissue_names, observed_tissues)
	return beta_vec, stderr_vec

# Print line to output file
def print_line_to_output_file_helper(t, snp_id, ensamble_id, pvalue_fe, pvalue_re2, i_square, pvalue_q, mvalue_vec, mvalue_stderr_vec, beta, stderr, num_observed_tissues):
	t.write(snp_id + ':' + ensamble_id + '\t' + num_observed_tissues + '\t' + pvalue_fe + '\t' + pvalue_re2 + '\t' + i_square + '\t' + pvalue_q)
	t.write('\t' + '\t'.join(mvalue_stderr_vec))
	t.write('\t' + '\t'.join(mvalue_vec))
	total_num_tissues = len(beta)
	for i in range(total_num_tissues):
		t.write('\t' + beta[i] + '\t' + stderr[i])
	t.write('\n')
	return t

eqtl_file = sys.argv[1]
organized_metatissue_output_dir = sys.argv[2]
metatissue_output_dir = sys.argv[3]
preprocess_metatissue_dir = sys.argv[4]
tissue_list_file = sys.argv[5]

# Extract vector of tissue names from tissue_list_file input file
tissue_names = extract_tissue_names(tissue_list_file)

# Open output file handle
output_file = organized_metatissue_output_dir + 'gemmatissue_output_v8.txt'
t = open(output_file, 'w')

# Print header line only to output file
t = print_header_to_output_file(t, tissue_names)

# Open file containing list of eqtls tested
f = open(eqtl_file)

# Loop through each variant gene pair tested
head_count = 0
for line in f:
	data = line.rstrip().split()
	if head_count == 0:  # Skip header
		head_count = head_count + 1
		continue
	# Extract variant and gene names
	snp_id = data[0]
	ensamble_id = data[1]

	# Extract dictionary list of all tissues observed for this variant gene pair from tissue info file
	tissue_info_file = preprocess_metatissue_dir + 'tissue_info_' + snp_id + '_' + ensamble_id + '.txt'
	observed_tissues = extract_observed_tissues_from_tissue_info_file(tissue_info_file)

	# Extract relevent fields from metasoft output file
	# This includes:
	#####1. PVALUE_FE
	#####2. PVALUE_RE2
	#####3. I_SQUARE
	#####4. PVALUE_Q
	#####5. vector of metatissue standard errors
	#####6. vector of metatissue mvalues
	metasoft_output_file = metatissue_output_dir + snp_id + '_' + ensamble_id + '_metatiss.metasoft.output.txt.gz'
	pvalue_fe, pvalue_re2, i_square, pvalue_q, mvalue_stderr_vec, mvalue_vec = extract_relevent_fields_from_metasoft_output(metasoft_output_file, tissue_names, observed_tissues)

	# Extract relevent fields from metatissue beta file
	# This includes
	#####1. Beta vector
	#####2. stderr vector
	metatissue_beta_file = metatissue_output_dir + snp_id + '_' + ensamble_id + '_metatiss.beta.std.txt.gz'
	beta, stderr = extract_relevent_fields_from_metatissue_beta_file(metatissue_beta_file, tissue_names, observed_tissues)

	# Print line to output file
	t = print_line_to_output_file_helper(t, snp_id, ensamble_id, pvalue_fe, pvalue_re2, i_square, pvalue_q, mvalue_vec, mvalue_stderr_vec, beta, stderr, str(len(observed_tissues)))

# Close filehandles
f.close()
t.close()
