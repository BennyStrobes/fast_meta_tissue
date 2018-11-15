import numpy as np 
import os
import sys
import random
import pdb



def get_tissues(file_name):
	f = open(file_name)
	arr = []
	for line in f:
		line = line.rstrip()
		arr.append(line)
	return np.asarray(arr)

# Open output file handles and write header
def open_handle(file_name):
	t = open(file_name,'w')
	t.write('snp_id\tensamble_id\ttop_tissue_name\ttissue_af\tpvalue\tfdr\n')
	return t


def get_gene_array_for_tissue(egene_file):
	f = open(egene_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(line)
	return arr

def get_autosomal_chroms():
	dicti = {}
	for chromer in range(1,23):
		dicti[str(chromer)] = 1
	return dicti

def get_egene(used_genes, egene_arr):
	autosomal_chroms = get_autosomal_chroms()
	converged = False
	while converged == False:
		liner = random.choice(egene_arr)
		gene_name = liner.split()[1]
		snp_id = liner.split()[0]
		snp_chrom_num = liner.split()[0].split('_')[0].split('hr')[1]
		#if gene_name not in used_genes and snp_chrom_num in autosomal_chroms and len(snp_id) <= 100:
		if gene_name not in used_genes:
			converged = True
			used_genes[gene_name] = 1
	return liner, used_genes


def randomly_select_egenes(cis_eqtl_file, tissues, trans_eqtl_file, cis_eqtl_bh_corrected_dir, cis_eqtl_egene_file_suffix):
	# Create gene array for each tissue
	tissue_specific_gene_array = {}
	for tissue in tissues:
		egene_file = cis_eqtl_bh_corrected_dir + tissue + cis_eqtl_egene_file_suffix
		tissue_specific_gene_array[tissue] = get_gene_array_for_tissue(egene_file)
	t = open_handle(cis_eqtl_file)
	f = open(trans_eqtl_file)
	used_genes = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[8]
		egene_line, used_genes = get_egene(used_genes, tissue_specific_gene_array[tissue])
		snp_id = egene_line.split()[0]
		ensamble_id = egene_line.split()[1]
		maf = egene_line.split()[3]
		pvalue = egene_line.split()[4]
		fdr = egene_line.split()[5]
		t.write(snp_id + '\t' + ensamble_id + '\t' + tissue_name + '\t' + maf + '\t' + pvalue + '\t' + fdr + '\n')
	f.close()
	t.close()
	return

trans_eqtl_file = sys.argv[1]
tissue_list_file = sys.argv[2]
cis_eqtl_bh_corrected_dir = sys.argv[3]
cis_eqtl_file = sys.argv[4]
top_cis_eqtl_file = sys.argv[5]





# Extract ordered arr of tissues
tissues = get_tissues(tissue_list_file)



randomly_select_egenes(cis_eqtl_file, tissues, trans_eqtl_file, cis_eqtl_bh_corrected_dir, '_egenes_and_random_variant_fdr_0.5.txt')


randomly_select_egenes(top_cis_eqtl_file, tissues, trans_eqtl_file, cis_eqtl_bh_corrected_dir, '_egenes_and_top_variant_fdr_0.5.txt')


# Loop through each tissue
#for tissue in tissues:
	# Cis eqtl file for this tissue
	#tissue_cis_eqtl_file = cis_eqtl_bh_corrected_dir + tissue + '_bh_corrected_0.5.txt'

	#top_egenes, random_egenes = extract_list_of_egenes(tissue_cis_eqtl_file)