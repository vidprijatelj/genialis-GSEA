#!/usr/bin/python3

import pandas as pd
import numpy as np
import os

#	set directory from which we work
cur_dir = os.getcwd()

#	Temp: set filename of our data
data_file_expressions = cur_dir + "/leukemia.txt"
data_pathways = cur_dir + "/pathways.txt"

#subset = np.genfromtxt(data_pathways, delimiter = "\t")
#print(subset)

#	Clean-up our pathways file as per Subramanian et al 2005
#	Any pathways/gene sets with less than 15 members are exluded,
#	as well as pathways/gene sets with "OBSOLETE" as one of the
#	qualifiers in their description



#	As we open a file we want to keep all data intact
#	hence header = None and index_col = None
df = pd.read_table(data_file_expressions, sep = '\t', header = None, index_col = None)

#	Define first row as a future header/column names
header = df.iloc[0,0:]

#	Remove first row entire
df = df[1:]

#	Define column names as our header
df.columns = header

#	Melt our dataframe into a manageable form
df = pd.melt(df, id_vars=['gene/patient'], var_name = 'phenotype')
#	print(df)

#	Parse through out dataframe and define a dictionary,
#	that holds our whole set of gene dictionaries
#	That dictionary -> {gene : {phenotype1 : [counts], 
#								phenotype2 : [counts]}}
def input_signals (input_dataframe):

	input_dict = dict()

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
				input_dict[gene][phenotype] = [count] # GOTTA DEFINE AS A LIST
		else:
	#	Define set[gene] as a new dictionary
	#	Define set[gene][phenotype] as a new dictionary
			input_dict[gene] = dict()
			input_dict[gene][phenotype] = dict()
			input_dict[gene][phenotype] = [count]

	return input_dict


input_dict = input_signals(df)

#	Function to calculate our signal to noise ratio
def signal_to_noise (list1, list2):
	mju1 = np.mean(list1)
	mju2 = np.mean(list2)
	sigma1 = np.std(list1)
	sigma2 = np.std(list2)
	return ( (mju1 - mju2) / (sigma1 + sigma2) )

#	list1 = set['TUBA1']['ALL']
#	list2 = set['TUBA1']['AML']
#	test = signal_to_noise(list1, list2)
#	print(test)

#	Function that parses through the dictionary of input signals
#	and returns dictionary in a form of gene:signal to noise ratio
def get_signal_to_noise (dictionary_in):
	dictionary_out = dict()
	for gene in dictionary_in:
		list1 = dictionary_in[gene]['ALL']
		list2 = dictionary_in[gene]['AML']
		dictionary_out[gene] = signal_to_noise(list1, list2)
	return dictionary_out

test = get_signal_to_noise(input_dict)

#	print(test)

#	Again we create a DataFramem, modify it and sort the values
def sort_values (dictionary_in):
	dataframe_out = pd.DataFrame(dictionary_in, index = [0])
	dataframe_out = dataframe_out.transpose()
	#	Define the values as Signal-to-Noise Ratio (StNR)
	dataframe_out.columns = ['StNR']
	dataframe_out = dataframe_out.sort_values(by = 'StNR', ascending = False)
	return dataframe_out

reference_dataframe = sort_values(test)
print(reference_dataframe)