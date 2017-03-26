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
#				output_file.write(line + '\n')
				gene_set_dict[split_line[0]] = split_line[2 : len(split_line)]
	return gene_set_dict


#	Clean-up our gene list file and melt the values
#	So we can parse through it
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
	return df



#	Function to scramble our phenotypes inside the dataframe
#	Work with beginning gene list!
#	Needed to create ES forward on
def scramble_gene_list_dataframe(input_dataframe):
	phenotype = input_dataframe['phenotype']
	phenotype = phenotype.sample(frac = 1).reset_index(drop = True)
	input_dataframe['phenotype'] = phenotype



#	Parse through out dataframe and define a new dataframe,
#	that holds our whole set of gene dictionaries
#	That dictionary -> {gene : {phenotype1 : [counts], 
#								phenotype2 : [counts]}}
def raw_gene_list_dataframe(input_dataframe):

	input_dict=dict()

	#	Iterate through rows of our input dataframe
	for row in input_dataframe.itertuples():
	#	set: {gene : {phenotype1 : count, phenotype2 : count} }
		gene = row[1]
		phenotype = row[2]
		count = int(row[3]) # EXPLICIT CALL AS INT!

	#	Check whether the key = gene already exist
		if gene in input_dict:
	#	Check whether the key = gene[phenotype] already exist
			if phenotype in input_dict[gene]:
	#	Append count value to the key gene[phenotype]
				input_dict[gene][phenotype].append(count)
			else:
	#	Define set[gene][phenotype] as a new dictionary
				input_dict[gene][phenotype] = dict()
				input_dict[gene][phenotype] = [count] # GOTTA DEFINE AS A LIST!
		else:
	#	Define set[gene] as a new dictionary
	#	Define set[gene][phenotype] as a new dictionary
			input_dict[gene] = dict()
			input_dict[gene][phenotype] = dict()
			input_dict[gene][phenotype] = [count]

	return input_dict



#	Function that parses through the dictionary of input signals
#	and returns dictionary in a form of gene:signal to noise ratio
#	in regards ALL on top and AML on bottom
def get_signal_to_noise(input_dictionary):
	dictionary_out = dict()

	#	Function to calculate our signal to noise ratio
	#	Called inside func get_signal_to_noise!
	def signal_to_noise(list1, list2):
		mju1 = np.mean(list1)
		mju2 = np.mean(list2)
		sigma1 = np.std(list1)
		sigma2 = np.std(list2)
		return ( (mju1 - mju2) / (sigma1 + sigma2) )


	for gene in input_dictionary:
		list1 = input_dictionary[gene]['ALL']
		list2 = input_dictionary[gene]['AML']
		dictionary_out[gene] = signal_to_noise(list1, list2)
	return dictionary_out



#	We take a dictionary, turn it into a dataframe, 
#	sort the dataframe in descending order in regards 
#	to of signal-to-noise ratio. High StNR is correlated
#	to ALL, low StNR is correlated to AML
	
def sort_gene_list_dataframe_values(input_dictionary):
	dataframe_out = pd.DataFrame(input_dictionary, index = [0])
	dataframe_out = dataframe_out.transpose()
	#	Define the values as Signal-to-Noise Ratio (StNR)
	dataframe_out.columns = ['StNR']
	dataframe_out = dataframe_out.sort_values(by = 'StNR', ascending = False)
	return dataframe_out




#	i = len(reference_dataframe.iloc[0:,0:0])
#	k = reference_dataframe.iloc[i-1:i, 0:1]
#	print(reference_dataframe.iloc[i-1:i, 0:0])
#	print(i, k)
#	print(reference_dataframe.index[0])
#	print(reference_dataframe.iloc[0]['StNR'])

#-----------------------------------------------------------------------------#




#gene_set_dict = gene_set_file_cleanup(data_pathways)
#beginning_gene_list_df = gene_list_dataframe_cleanup(data_gene_expressions)
#working_gene_list_df = raw_gene_list_dataframe(beginning_gene_list_df)
#working_gene_list_dict = get_signal_to_noise(working_gene_list_df)
#working_gene_list_df = sort_gene_list_dataframe_values(working_gene_list_dict)

#scramble_gene_list_dataframe(beginning_gene_list_df)
#print(beginning_gene_list_df)
#print(working_gene_list_df)



#------------------------------------------------------------------------------#
#	Calculate ES by creating Pandas tables
#------------------------------------------------------------------------------#

def create_ES_tables(file_name_expressions, file_name_pathways):

	gene_set_dict = gene_set_file_cleanup(file_name_pathways)
	beginning_gene_list_df = gene_list_dataframe_cleanup(file_name_expressions)
	working_gene_list_df = raw_gene_list_dataframe(beginning_gene_list_df)
	working_gene_list_dict = get_signal_to_noise(working_gene_list_df)
	working_gene_list_df = sort_gene_list_dataframe_values(working_gene_list_dict)

	#	define N
	N = working_gene_list_df.shape[0]
	test_dataframe = pd.DataFrame()
#	working_gene_list.df = working_gene_list_df.round(6)
#	print(working_gene_list.df)

	print('\nCalculating ES\n')

	#	Define list gene_set_series that will serve as an index in our
	#	final dataframe
	#	Define dictionary ES_dictionary that will serve as columns for
	#	all our ES's (original and permutated ones) in our final dataframe

#	with open(file_name_pathways+'_clean.txt', 'r') as gene_set:
#		for line in gene_set:
#			split_line = re.split(r'\t+', line)
#			pathway = split_line[0]
#			gene_set_series.append(pathway)

	#	Get our ES_series and append it to our dataframe

	#	1000 permutations in order to create ESnull	
	for i in range(0, 10):	

		ESlist = list()
		ESdict = dict()

		for pathway in gene_set_dict:

			#print(pathway, gene_set_dict[pathway])
			#	Parse through whole data in order to calculate ES for
			#	each gene set

				#	Get ES for each gene set
			def get_ES(pathway):

				#	Calculate Nr and define Nh
				#	Iterating through row by row				
				Nh = 0
				Nr = 0
				for row in working_gene_list_df.itertuples():
					gene = str(row[0])
					r = abs(row[1])
					if gene in gene_set_dict[pathway]:
						Nr = Nr + r
						Nh = Nh + 1
					
				#Define Pmiss
				Pmiss = 1 / (N - Nh)

				ES = 0
				ESlist = list()
				#	Get all ES values
				for row in working_gene_list_df.itertuples():
					gene = str(row[0])
					r = abs(row[1])
					if gene in gene_set_dict[pathway]:
						ES = ES + (r / Nr)
					else:
						ES = ES - Pmiss
					ESlist.append(ES)

				if abs(max(ESlist)) > abs(min(ESlist)):
					ES = max(ESlist)
				else:
					ES = min(ESlist)
				return(ES)

			ES = get_ES(pathway)
			ESdict[pathway] = ES

		#	Create a permutation and start over again

		index = list(ESdict.keys())
		values = list(ESdict.values())

		test_series = pd.Series.to_frame(pd.Series(ESdict, name='ES'+str(i)))
		if i == 0:
			test_dataframe = test_series
		else:
			test_dataframe = pd.merge(test_dataframe, test_series, left_index = True, right_index = True)
		#else:
		#	test_dataframe.append(test_series)

		scramble_gene_list_dataframe(beginning_gene_list_df)
		working_gene_list_df = raw_gene_list_dataframe(beginning_gene_list_df)
		working_gene_list_dict = get_signal_to_noise(working_gene_list_df)
		working_gene_list_df = sort_gene_list_dataframe_values(working_gene_list_dict)

		
#		print(gene_set_series)
#		print(ES_series)
	print(test_dataframe)



#------------------------------------------------------------------------------#
#	Calculate ES by using dictionaries
#------------------------------------------------------------------------------#

def create_ES_dictionaries(file_name_expressions, file_name_pathways):

	gene_set_dict = gene_set_file_cleanup(file_name_pathways)
	beginning_gene_list_df = gene_list_dataframe_cleanup(file_name_expressions)
	working_gene_list_df = raw_gene_list_dataframe(beginning_gene_list_df)
	working_gene_list_dict = get_signal_to_noise(working_gene_list_df)
	working_gene_list_df = sort_gene_list_dataframe_values(working_gene_list_dict)

	#	define N
	N = working_gene_list_df.shape[0]
	test_dictionary = dict()
#	working_gene_list.df = working_gene_list_df.round(6)
#	print(working_gene_list.df)

	print('\nCalculating ES\n')

	#	Define list gene_set_series that will serve as an index in our
	#	final dataframe
	#	Define dictionary ES_dictionary that will serve as columns for
	#	all our ES's (original and permutated ones) in our final dataframe

#	with open(file_name_pathways+'_clean.txt', 'r') as gene_set:
#		for line in gene_set:
#			split_line = re.split(r'\t+', line)
#			pathway = split_line[0]
#			gene_set_series.append(pathway)

	#	Get our ES_series and append it to our dataframe

	#	1000 permutations in order to create ESnull	
	for pathway in gene_set_dict:
		test_dictionary[pathway] = dict()

	for i in range(0, 10):	

		ESlist = list()

		for pathway in gene_set_dict:

			#print(pathway, gene_set_dict[pathway])
			#	Parse through whole data in order to calculate ES for
			#	each gene set

				#	Get ES for each gene set
			def get_ES(pathway):

				#	Calculate Nr and define Nh
				#	Iterating through row by row				
				Nh = 0
				Nr = 0
				for row in working_gene_list_df.itertuples():
					gene = str(row[0])
					r = abs(row[1])
					if gene in gene_set_dict[pathway]:
						Nr = Nr + r
						Nh = Nh + 1
					
				#Define Pmiss
				Pmiss = 1 / (N - Nh)

				ES = 0
				ESlist = list()
				#	Get all ES values
				for row in working_gene_list_df.itertuples():
					gene = str(row[0])
					r = abs(row[1])
					if gene in gene_set_dict[pathway]:
						ES = ES + (r / Nr)
					else:
						ES = ES - Pmiss
					ESlist.append(ES)

				if abs(max(ESlist)) > abs(min(ESlist)):
					ES = max(ESlist)
				else:
					ES = min(ESlist)
				return(ES)

			ES = get_ES(pathway)
			test_dictionary[pathway]['ES'+str(i)] = ES

		#	Create a permutation and start over again

		scramble_gene_list_dataframe(beginning_gene_list_df)
		working_gene_list_df = raw_gene_list_dataframe(beginning_gene_list_df)
		working_gene_list_dict = get_signal_to_noise(working_gene_list_df)
		working_gene_list_df = sort_gene_list_dataframe_values(working_gene_list_dict)

		
#		print(gene_set_series)
#		print(ES_series)
	print(test_dictionary)

t = time.process_time()
create_ES_dictionaries(data_gene_expressions, data_pathways)
elapsed_time = time.process_time() - t
print(elapsed_time)