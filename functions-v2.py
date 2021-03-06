#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import re
import time



#	set directory from which we work
cur_dir = os.getcwd()

#	Temp: set filename of our data
data_gene_expressions = cur_dir + "/leukemia.txt"
data_pathways = cur_dir + "/pathways.txt"



#	Clean-up our gene sets file as per Subramanian et al 2005.
#	Any pathways/gene sets with less than 15 genes are exluded,
#	as well as pathways/gene sets with "OBSOLETE" as one of the
#	qualifiers in their description.
#	We output the gene set in a new dict
def gene_set_file_cleanup(file_name):
	gene_set_dict = dict()
	with open(file_name, 'r') as temp_file, open(file_name+'_clean.txt', 'w') as output_file:
		#firstNlines=temp_file.readlines()[0:5]

		print('\nFollowing gene sets will be removed: \n')
		for line in temp_file:
			#	Remove \n from end of line
			line = line.rstrip('\n')
			split_line = re.split(r'\t+', line)			
			obsolete = "OBSOLETE"
			if ( (obsolete in split_line[1]) == True ) or ( len(split_line) - 2 < 15 ) :
				print('Removed:', split_line[0], '\t', 'value:', len(split_line) - 2)
			else:
				output_file.write(line + '\n')
				gene_set_dict[split_line[0]] = split_line[2 : len(split_line)]
	return gene_set_dict


#######################################################################

#--------#
#	Version 1
#	Dataframe methods only
#	~3s per permutation
#--------#
'''

#	Clean-up our gene list file and melt the values
#	So we can parse through it
def gene_list_dataframe_cleanup(file_name):
	#	As we open a file we want to keep all data intact
	#	hence header = None and index_col = None
	df = pd.read_table(file_name, sep = '\t', header = None, index_col = None)

	#	Define first row as a future header/column names
	header = df.iloc[0,1:]
	index = df.iloc[1:,0]
	#	Remove first row entire
	df = df[1:]
	df = df.drop(0, axis = 1)
	#	Define column names as our header
	df.columns = header
	df = df.set_index(index)
	df.columns.name = 'gene/patient'
	

	#	Melt our dataframe into a manageable form
	#df = pd.melt(df, id_vars=['gene/patient'], var_name = 'phenotype')
	return df



#	Function to scramble our phenotypes inside the dataframe.
#	Works with beginning gene list! Needed to create ES forward on
def scramble_gene_list_dataframe(input_dataframe):
	columns = list(input_dataframe)
	np.random.shuffle(columns)
	input_dataframe.columns = columns
	return input_dataframe



#	Parse through our dataframe and define a new dataframe,
#	where index = gene and columns = Signal-to-Noise Ratio
#	Returns a new dataframe
def signal_to_noise_sort_dataframe(input_dataframe):
	output_dict = dict()

	#	Function to calculate our signal to noise ratio	
	def signal_to_noise(list1, list2):
		mju1 = np.mean(list1)
		mju2 = np.mean(list2)
		sigma1 = np.std(list1)
		sigma2 = np.std(list2)
		StNR = ((mju1 - mju2) / (sigma1 + sigma2))
		return StNR
	
	#	Define unique genes in our dataframe
	unique_genes = input_dataframe.index.values

	#	Zip through our list of unique genes
	for gene in list(zip(unique_genes)):
		gene = gene[0]
		t = input_dataframe.loc[gene]
		ALL = t.loc['ALL'].values.tolist()
		ALL = list(map(int, ALL))
		AML = t.loc['AML'].values.tolist()
		AML = list(map(int, AML))
		StNR = signal_to_noise(ALL, AML)
		output_dict[gene] = StNR

	output_dataframe = pd.DataFrame.from_dict(output_dict, orient = 'index')
	output_dataframe.columns = ['StNR']
	output_dataframe = output_dataframe.sort_values(by = 'StNR', ascending = False)
	return output_dataframe
'''

#######################################################################


#--------#
#	Version 2
#	Mixed bag
#	~1.5s per permutation
#--------#

def gene_list_dataframe_cleanup(file_name):
	#	As we open a file we want to keep all data intact
	#	hence header = None and index_col = None
	df = pd.read_table(file_name, sep = '\t', header = None, index_col = None)

	#	Define first row as a future header/column names
	header = df.iloc[0,0:]

	#	Remove first row entire
	df = df[1:]

	#	Define column names as our header
	df.columns = header

	#	Melt our dataframe into a manageable form
	df = pd.melt(df, id_vars=['gene/patient'], var_name = 'phenotype')
	df['phenotype'] = df['phenotype'].str.replace('\n', '')
	return df



#	Function to scramble our phenotypes inside the dataframe
#	Work with beginning gene list!
#	Needed to create ES forward on
def scramble_gene_list_dataframe(input_dataframe):
	phenotype = input_dataframe['phenotype']
	phenotype = phenotype.sample(frac = 1).reset_index(drop = True)
	input_dataframe['phenotype'] = phenotype



#	Parse through out dataframe and define a new dictionary,
#	that holds our whole set of gene dictionaries
#	That dictionary -> {gene : {phenotype1 : [values], 
#								phenotype2 : [values]}}
#	Then we get the signal-to-noise ratio for each of our genes
#	and create an output sorted dataframe (sort by StNR descending)
def signal_to_noise_sort_dataframe(input_dataframe):

	input_dict = dict()
	output_dict = dict()

	#	Iterate through rows of our input dataframe
	for zip_tuple in list(zip(input_dataframe['gene/patient'], input_dataframe['phenotype'], input_dataframe['value'])):
	#	set: {gene : {phenotype1 : value, phenotype2 : value} }
		gene = zip_tuple[0]
		phenotype = zip_tuple[1].rstrip('\n') # HAVE TO STRIP \n!
		value = int(zip_tuple[2]) # EXPLICIT CALL AS INT!

	#	Check whether the key = gene already exist
		if gene in input_dict:
	#	Check whether the key = gene[phenotype] already exist
			if phenotype in input_dict[gene]:
	#	Append value value to the key gene[phenotype]
				input_dict[gene][phenotype].append(value)
			else:
	#	Define set[gene][phenotype] as a new dictionary
				input_dict[gene][phenotype] = dict()
				input_dict[gene][phenotype] = [value] # GOTTA DEFINE AS A LIST!
		else:
	#	Define set[gene] as a new dictionary
	#	Define set[gene][phenotype] as a new dictionary
			input_dict[gene] = dict()
			input_dict[gene][phenotype] = dict()
			input_dict[gene][phenotype] = [value]

	#	Function to calculate our signal to noise ratio
	def signal_to_noise(list1, list2):
		mju1 = np.mean(list1)
		mju2 = np.mean(list2)
		sigma1 = np.std(list1)
		sigma2 = np.std(list2)
		return ( (mju1 - mju2) / (sigma1 + sigma2) )	

	#	Create a new dictionary with key = gene : value = StNR
	for gene in input_dict:
		list1 = input_dict[gene]['ALL']
		list2 = input_dict[gene]['AML']
		output_dict[gene] = signal_to_noise(list1, list2)

	#	Return an output dataframe
	dataframe_out = pd.DataFrame.from_dict(output_dict, orient = 'index')
	dataframe_out.columns = ['StNR']
	dataframe_out = dataframe_out.sort_values(by = 'StNR', ascending = False)

	return dataframe_out



#	Get ES for each gene set in a pathway
def get_ES(input_dataframe, N, gene_set_dict, pathway):

	#	Create a copy of a working df with indeces that
	#	appear only in gene set of a current pathway.
	#	If a gene is not in a working df, pandas will 
	#	value it as a "NaN", hence we need to drop it
	present_genes_df = input_dataframe.copy()
	present_genes_df = present_genes_df.loc[gene_set_dict[pathway]].dropna(axis = 0)
	#	r = abs!
	present_genes_df = present_genes_df.apply(abs)

	#	Define Nh and Nr
	Nh = present_genes_df.shape[0]
	Nr = sum(present_genes_df.values)

	#	Positive values we add are always r / Nr!
	present_genes_df = present_genes_df.div(Nr)
		
	#Define Pmiss
	Pmiss = -(1 / (N - Nh))

	#	Create a copy of a working df with all values
	#	defined as Pmiss. We then merge both newly
	#	created df's.
	temporary_df = input_dataframe.copy()
	temporary_df[:] = Pmiss
	temporary_df.update(present_genes_df)

	ES = 0
	ES_list = list()

	#	Get all ES values
	#	We use the numpy cumsum method
	ES_list = np.cumsum(temporary_df['StNR'])
				
	#	Copy ES_list and make all values absolute
	ES_list_abs = (np.absolute(ES_list))
	#	Find an index of the biggest value in list of ES_list absolutes
	index_max_dev = np.argmax(ES_list_abs)
	#	ES is the maximum absolute value in a list of ES_list absolutes
	#	That way we can ignore pos and neg values.
	ES = ES_list[index_max_dev]
	return ES



#------------------------------------------------------------------------------#
#	Calculate ES by creating Pandas tables
#------------------------------------------------------------------------------#

def create_ES_tables(file_name_expressions, file_name_pathways):

	gene_set_dict = gene_set_file_cleanup(file_name_pathways)
	beginning_gene_list_df = gene_list_dataframe_cleanup(file_name_expressions).copy()
	working_gene_list_df = signal_to_noise_sort_dataframe(beginning_gene_list_df)

	ES_dictionary = dict()
	ES_dataframe = pd.DataFrame()

	#	define N
	N = working_gene_list_df.shape[0]


	print('\nCalculating ES\n')

	#	1000 permutations in order to create ESnull	
	for i in range(0, 200):	

		print('Creating permutation No.', + i)
		#	Parse through whole pathways data in order to calculate ES for
		#	each gene set (each pathway)
		for pathway in gene_set_dict:

			ES = get_ES(working_gene_list_df, N, gene_set_dict, pathway)
			ES_dictionary[pathway] = ES

		index = list(ES_dictionary.keys())
		values = list(ES_dictionary.values())

		current_ES_dataframe = pd.DataFrame.from_dict(ES_dictionary, orient = 'index')
		current_ES_dataframe.columns = ['ES' + str(i)]
		ES_dataframe = pd.concat([ES_dataframe, current_ES_dataframe], axis = 1)
	

		#	Create a permutation and start over again
		scramble_gene_list_dataframe(beginning_gene_list_df)
		working_gene_list_df = signal_to_noise_sort_dataframe(beginning_gene_list_df)

		
#		print(gene_set_series)
#		print(ES_series)
	print(ES_dataframe)



t = time.process_time()
create_ES_tables(data_gene_expressions, data_pathways)
elapsed_time = time.process_time() - t
print(elapsed_time)