import numpy as np 
import sys
import os
import pdb








# Extract dictionary list of all cis eqtls at FDR <= fdr_thresh
def extract_cis_eqtls(cis_eqtl_file):
	cis_eqtls = {}
	f = open(cis_eqtl_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			header = line
			continue
		snp_id = data[0]
		ensamble_id = data[1]
		cis_eqtls[snp_id] = line
	f.close()
	return cis_eqtls, header






tissue_name = sys.argv[1]
cis_eqtl_bh_corrected_dir = sys.argv[2]
fdr_thresh = float(sys.argv[3])
trans_eqtl_file = sys.argv[4]


# Extract dictionary list of all cis eqtls at FDR <= fdr_thresh
cis_eqtl_file = cis_eqtl_bh_corrected_dir + tissue_name + '_bh_corrected_' + str(fdr_thresh) + '.txt'
cis_eqtls, header = extract_cis_eqtls(cis_eqtl_file)

############
# Loop through trans eqtls. If variant overlaps, print line
output_file = cis_eqtl_bh_corrected_dir + tissue_name + '_cis_trans_overlaps_cis_eqtls_' + str(fdr_thresh) + '.txt'
t = open(output_file, 'w')
t.write(header + '\n')
f = open(trans_eqtl_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	trans_eqtl_fdr = float(data[-1])
	trans_eqtl_tissue = data[-4]
	if trans_eqtl_fdr > fdr_thresh or trans_eqtl_tissue != tissue_name:
		continue
	snp_id = data[0]
	if snp_id in cis_eqtls:
		t.write(cis_eqtls[snp_id] + '\n')
t.close()
