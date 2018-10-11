#Benjamin Strober
#bstrober3@gmail.com
#1/12/17
#A faster implementation of MetaTissue (Su et al) via utilization of GEMMA to estimate variance components as opposed EMMA.
import numpy as np
import os
import sys
import gzip
import pdb
from sklearn import linear_model


##################################################################################################################################
#run_fast_meta_tissue is the main/driver function for this module.
#This function is intended for 1 snp-gene pair at a time. So every variant will need its own genotype file. And every gene will need its own expression file.
#This function assumes all covariates are regressed out. See note in README concerning regressing out covariates in MetaTissue.
#The parameters are as follows:
###raw_output_dir --> Absolute directory to save results in.
###gene_list --> Text file where each line is the location of a specific tissues expression file (see http://genetics.cs.ucla.edu/metatissue/install_step2.html 'gene expression list file')
###genotype_file --> Tab deliminated file containing genotype (where non-tissue-specific covariates are regressed out) for all samples (see http://genetics.cs.ucla.edu/metatissue/install_step2.html EIGENSTRAT genotype file)
###tissue_info_file --> File containing information concerning which samples are expressed in which tissues for this particular gene. (see http://genetics.cs.ucla.edu/metatissue/install_step2.html tissue information file)
###heuristic --> Boolean variable (True or False) dictating whether to run the heuristic option. Should be set to False for now
###snp --> string of name SNPID
###gene --> string of geneID
###gemma_directory --> Absolute directory that contains gemma executable. I've been using GEMMA_95_alpha. I have included it in the 'packages' directory
###java_compiler --> Location of java compiler
###metasoft_jar --> location of metasoft_jar file. This is found in the downloaded metatissue package. I have also included it in the 'packages' directory
###han_eskin_pvalue_table --> location of this pvalue table. Found in the downloaded metatissue package. I have also included it in the 'packages' directory
###mt_pvalue_table --> location of this pvalue table. Found in the dowloaded metatissue package.  I have also included it in the 'packages' directory.
#Return a binary variable. True if more than one tissue express this snp-gene pair. False if only one tissue expresses this snp gene pair (not enough tissues to run meta-analysis on, duhhhhhh)
##################################################################################################################################
def run_fast_meta_tissue(raw_output_dir,gene_list_file,genotype_file,tissue_info_file,heuristic,snp,gene,gemma_directory,java_compiler,metasoft_jar,han_eskin_pvalue_table,mt_pvalue_table):
    output_dir = raw_output_dir + snp + '_' + gene + '_'
    #preprocess data to get in appropriate format for GEMMA 
    run_binary,model_data = preprocess_fast_meta_tissue(output_dir,gene_list_file,genotype_file,tissue_info_file)
    #If only one tissue, break out of this function
    if run_binary == False:
        return False
    #The heuristic ratio for each tissue. See MetaTissue suppliment. If heuristic == False, this will return vector of ones (ie. ratio has no effect)
    heuristic_ratio = heuristic_option(output_dir,model_data,heuristic)
    #Use GEMMA to extract variance components of GLMM
    estimate_variance_component_command = gemma_directory + 'gemma -p ' + output_dir + 'phenotype.txt -k ' + output_dir + 'relatedness_matrix.txt -vc 1 -outdir ' + raw_output_dir + ' -o ' + snp + '_' + gene + '_' + 'gemma_vc_output'
    os.system(estimate_variance_component_command)
    #Fit the model
    fit_model(output_dir, model_data,heuristic_ratio,snp,gene,heuristic,java_compiler,metasoft_jar,han_eskin_pvalue_table,mt_pvalue_table)
    return True



def create_relatedness_matrix(sample_vec,n_measurements,output_file):
    temp = np.zeros((n_measurements,n_measurements)).astype(str)
    for i,val in enumerate(sample_vec):
        indices = np.where(sample_vec == val)[0]
        for index in indices:
            temp[index,i] = '1.0'
            temp[i,index] = '1.0'
    print_matrix_line_by_line(output_file,temp,n_measurements)
    return temp
def print_matrix_line_by_line(output_file,temp,n):
    t = open(output_file,'w')
    for i in range(n):
        t.write('\t'.join(temp[i,:]) + '\n')
    t.close()

def extract_valid_indices(data, missing_donors):
    indices = []
    for i, val in enumerate(data):
        if val not in missing_donors:
            indices.append(i)
    return indices

def make_phenotype_vector(gene_list_file,output_file, missing_donors):
    phen_vec = np.asarray([])
    sample_vec = np.asarray([])
    tissue_counts = []
    f = open(gene_list_file)
    for itera,line in enumerate(f):
        filer = line.rstrip()
        g = open(filer)
        count = 0
        for line in g:
            line = line.rstrip()
            data = np.asarray(line.split())
            if count == 0:
                count = count + 1
                valid_indices = extract_valid_indices(data, missing_donors)
                sample_vec = np.hstack((sample_vec,data[valid_indices]))
                continue
            phen_vec = np.hstack((phen_vec,data[valid_indices]))
            tissue_counts.append(len(data[valid_indices]))
    np.savetxt(output_file,phen_vec,delimiter='\n',fmt="%s")
    return tissue_counts,sample_vec,phen_vec
def print_design_matrix_helper(output_file1,output_file2,design_matrix):
    t = open(output_file1,'w')
    row,col = design_matrix.shape 
    temp = design_matrix[:,0:-1]
    for i in range(row):
        t.write('\t'.join(temp[i,:]) + '\n')
    t.close()
    t = open(output_file2,'w')
    t.write('rs1, A, T, ' + ', '.join(design_matrix[:,-1]))
    t.close()



#The heuristic ratio for each tissue. See MetaTissue suppliment. If heuristic == False, this will return vector of ones (ie. ratio has no effect)
def heuristic_option(output_dir,model_data,heuristic):
    if heuristic == True:
        combined_std = get_combined_std_heuristic(output_dir,model_data)
        tbt_std = get_tbt_std_heuristic(output_dir,model_data)
        heuristic_ratio = np.asarray(tbt_std/combined_std)
    else:
        heuristic_ratio = np.ones(len(model_data[3]))
    return heuristic_ratio

def get_combined_std_heuristic(output_dir,model_data):
    y = model_data[0]
    x = model_data[1]
    model = linear_model.LinearRegression(fit_intercept=False)
    modelfit = model.fit(x,y)
    mse = np.mean(np.square((y-modelfit.predict(x))))
    var = np.diagonal(np.linalg.inv(np.dot(np.transpose(x),x))*mse)
    return np.sqrt(var)[len(var)/2:]
def get_tbt_std_heuristic(output_dir,model_data):
    y = model_data[0]
    x = model_data[1]
    tissue_counts = model_data[3]
    sdev = []
    cumulative=0
    for tiss_index,tiss_count in enumerate(tissue_counts):
        temp_y = y[cumulative:cumulative+tiss_count]
        temp_x = np.hstack((x[cumulative:cumulative+tiss_count,tiss_index:tiss_index+1],x[cumulative:cumulative+tiss_count,len(tissue_counts) + tiss_index:len(tissue_counts) + tiss_index + 1]))
        model = linear_model.LinearRegression(fit_intercept=False)
        modelfit = model.fit(temp_x,temp_y)
        mse = np.mean(np.square((temp_y-modelfit.predict(temp_x))))
        var = np.diagonal(np.linalg.inv(np.dot(np.transpose(temp_x),temp_x))*mse)   
        sdev.append(np.sqrt(var)[1])
        cumulative=cumulative+tiss_count
    return sdev

# First need to remove rows of tissue_info_file (ie. donors) that don't have observed genotype data
def remove_missing_donors(temp,missing_donors):
    good_rows = []
    for row_num in range(temp.shape[0]):
        if temp[row_num,0] not in missing_donors:
            good_rows.append(row_num)
    temp = temp[good_rows,:]
    return temp

def find_samples_in_each_tissue(input_file, missing_donors):
    temp = np.loadtxt(input_file,delimiter='\t',dtype=str,comments=' ')
    # First need to remove rows of tissue_info_file (ie. donors) that don't have observed genotype data
    temp = remove_missing_donors(temp,missing_donors)
    row,col = temp.shape 
    n_tiss = col -1
    n_samp = row-1
    index_to_tissue = temp[0,1:]
    index_to_array = {}
    for index in range(n_tiss):
        index_to_array[index] = np.where(temp[1:,1+index] == '1')[0]
    return n_samp,n_tiss,index_to_tissue,index_to_array

def remove_nas_from_vector(geno):
    new_geno = []
    for ele in geno:
        if ele == 'NA':
            continue
        new_geno.append(ele)
    return(np.asarray(new_geno))

#Create design matrix
def make_design_matrix(tissue_counts,genotype_file,n_tiss,index_to_array,n_measurements):
    mat = np.zeros((n_measurements,2*n_tiss)).astype(str)
    cumulative = 0
    geno = np.loadtxt(genotype_file,delimiter='\t',dtype=str)
    geno = remove_nas_from_vector(geno)
    for itera,count in enumerate(tissue_counts):
        mat[cumulative:(cumulative+count),itera] = np.ones(count)
        indices = index_to_array[itera]
        mat[cumulative:(cumulative+count),itera+n_tiss] = geno[indices]
        cumulative = cumulative + count 
    return mat

# Get list of donors that are missing genotype data for this variant
def extract_donors_with_missing_genotype_data(genotype_file, tissue_info_file):
    # First create list of the ordered names of all the donors
    all_donors = []
    f = open(tissue_info_file)
    head_count = 0 # Skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0: # Skip header
            head_count = head_count + 1
            continue
        # Add donor name to array
        all_donors.append(data[0])
    f.close()

    # Using the ordered list of all the donors, go into the genotype file to see which donors are NA
    missing_donors = {}
    f = open(genotype_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        genotype_vec = data
    f.close()
    # Check to make sure len(all_donors) == len(genotype_vec)
    if len(all_donors) != len(genotype_vec):
        print('fundamental assumption error!')
    for i, dosage in enumerate(genotype_vec):
        if dosage == 'NA':
            missing_donors[all_donors[i]] = 1
    return missing_donors


##Driver to preprocess/reformat data into correct format for GEMMA
def preprocess_fast_meta_tissue(output_dir,gene_list_file,genotype_file,tissue_info_file):
    print('Reformatting data...')
    model_data = []

    ##########################################
    # Get list of donors that are missing genotype data for this variant
    missing_donors = extract_donors_with_missing_genotype_data(genotype_file, tissue_info_file)

    #############################################
    #Find exactly which samples are expressed in which tissue
    n_samp,n_tiss,index_to_tissue,index_to_array = find_samples_in_each_tissue(tissue_info_file, missing_donors)

    #If number of tissues expressed == 1, get out.
    if n_tiss == 1:
        return False,model_data
    #Make phenotype vector, y
    tissue_counts,sample_vec,y = make_phenotype_vector(gene_list_file,output_dir + 'phenotype.txt', missing_donors)
    n_measurements = len(sample_vec)
    ################################
    #Save the relatedness matrix. This is a big file and takes a while to save
    #If anyone knows how to speed this up, it'd definitely help
    d = create_relatedness_matrix(sample_vec,n_measurements,output_dir + 'relatedness_matrix.txt')  
    ###################################
    design_matrix = make_design_matrix(tissue_counts,genotype_file,n_tiss,index_to_array,n_measurements)
    print_design_matrix_helper(output_dir + 'design_matrix_1.txt',output_dir + 'design_matrix_2.bim',design_matrix)
    #Return a data vector composed of:
    ##phenotype vector
    ##Design matrix
    ##Relateness matrix
    ##Number of samples in each tissue
    model_data = [np.transpose(np.asmatrix(y.astype(float))),design_matrix.astype(float),d.astype(float),tissue_counts]
    print('Finished reformatting data: ')
    return True,model_data
#Parse GEMMA output file and return variance component estimates
def extract_variance_coef_from_output_text(file_name):
    f = open(file_name)
    binary = False
    for line in f:
        if line.startswith('## sigma2 estimates'):
            line = line.rstrip()
            data = line.split()
            var1 = float(data[-2])
            var2 = float(data[-1])
            binary = True
    if binary == False:
        print('Error in GEMMA OUTPUT')
    return var1,var2

#Take a covariance matrix and convert it into a correlation matrix
def convert_covariance_to_correlation(mat):
    row,col = mat.shape
    new_mat = np.identity(row)
    for i in range(row):
        for j in range((i+1),row):
            corr = mat[i,j]/np.sqrt(mat[i,i]*mat[j,j])
            new_mat[i,j] = corr
            new_mat[j,i] = corr
    return new_mat
def print_helper_fast_meta_beta_std(output_file,snp,gene,beta,std):
    t = open(output_file,'w')
    t.write(snp + ':' + gene)
    for i,val in enumerate(beta):
        t.write(' ' + str(val[0,0]) + ' ' + str(std[i]))
    t.close()


def print_helper_fast_meta_corr(output_file,corr):
    t = open(output_file,'w')
    row,col = corr.shape
    for j in range(col):
        for i in range(row):
            if i == 0 and j == 0:
                t.write(str(corr[i,j]))
            else:
                t.write(' ' + str(corr[i,j]))
    t.close()

def print_helper_fast_meta_sigmag(output_file,var1,var2):
    fraction = var1/(var1+var2)
    t = open(output_file,'w')
    t.write(str(fraction))
    t.close()

def print_helper_fast_meta_metasoft_script(output_file,output_dir,tissue_counts,heuristic,java_compiler,metasoft_jar,han_eskin_pvalue_table,mt_pvalue_table):
    t = open(output_file,'w')
    t.write('#/bin/bash\n')
    t.write(java_compiler + ' -Xmx1G -jar ' + metasoft_jar + ' -pvalue_table ' + han_eskin_pvalue_table + ' -input ')
    t.write(output_dir + 'metatiss.beta.std.txt.gz -output ' + output_dir + 'metatiss.metasoft.output.txt.gz -correlation ')
    t.write(output_dir + 'metatiss.corr.txt.gz -log ' + output_dir + 'metatiss.metasoft.log.txt')
    if heuristic == True:
        stringer = ' '.join(np.asarray(tissue_counts).astype(str))
        t.write(' -sample_per_study "' + stringer + '"' + ' -mt_pvalue_table ' + mt_pvalue_table + ' -sigmag ' + output_dir + 'metatiss.sigmag.txt.gz')
    t.write(' -mvalue -mvalue_prior_sigma 0.4 -mvalue_p_thres 1.0 -seed 1 -mvalue_method mcmc')
    t.close()

#Print output files from glmm into correct format for input to meta_soft
def print_helper_fast_meta(output_dir,snp,gene,beta,updated_sdev,corr,var1,var2,model_data,heuristic,java_compiler,metasoft_jar,han_eskin_pvalue_table,mt_pvalue_table):
    print_helper_fast_meta_beta_std(output_dir + 'metatiss.beta.std.txt',snp,gene,beta,updated_sdev)
    print_helper_fast_meta_corr(output_dir + 'metatiss.corr.txt',corr)
    print_helper_fast_meta_sigmag(output_dir + 'metatiss.sigmag.txt',var1,var2)
    print_helper_fast_meta_metasoft_script(output_dir + 'meta_soft.sh',output_dir,model_data[3],heuristic,java_compiler,metasoft_jar,han_eskin_pvalue_table,mt_pvalue_table)
    os.system('gzip ' + output_dir + 'metatiss.beta.std.txt')
    os.system('gzip ' + output_dir + 'metatiss.corr.txt')
    os.system('gzip ' + output_dir + 'metatiss.sigmag.txt')

#Run metasoft
def run_meta_soft(output_dir):
    meta_soft_command = 'sh ' + output_dir + 'meta_soft.sh'
    os.system(meta_soft_command)

def fit_model(output_dir,model_data,heuristic_ratio,snp,gene,heuristic,java_compiler,metasoft_jar,han_eskin_pvalue_table,mt_pvalue_table):
    y = model_data[0] # phenotype vector
    x = model_data[1] # design matrix
    dd = model_data[2] #relatedness_matrix
    #parse gemma output file to get variance component estimates
    print('extracting')
    var1,var2 = extract_variance_coef_from_output_text(output_dir + 'gemma_vc_output.log.txt')
    print('create linear model')
    #Compute linear model
    Sigma = var1*dd + np.identity(len(y))*var2
    print('taking inverse')
    inv_sig = np.linalg.inv(Sigma)
    print('taking inverse pt 2')
    variances = np.linalg.inv(np.dot(np.dot(np.transpose(x),inv_sig),x)) #Covariance matrix of learned betas
    print('done taking inv')
    non_intercept_length = variances.shape[0]/2 #Number of tissues
    beta = np.dot(np.dot(np.dot(variances,np.transpose(x)),inv_sig),y)[non_intercept_length:] #learned betas
    sdev = np.sqrt(np.diag(variances))[non_intercept_length:] #learned standard deviations
    corr = convert_covariance_to_correlation(variances)[non_intercept_length:,non_intercept_length:] #correlation matrix
    updated_sdev = sdev*heuristic_ratio #update sdev with heuristic ratio. Recall if heuristic == False, stdev will stay the same
    print('printing')
    #print output files/model into correct format for input to meta_soft
    print_helper_fast_meta(output_dir,snp,gene,beta,updated_sdev,corr,var1,var2,model_data,heuristic,java_compiler,metasoft_jar,han_eskin_pvalue_table,mt_pvalue_table) 
    #Run meta_soft on output from model
    run_meta_soft(output_dir)


