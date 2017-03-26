#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import re
import time

np.set_printoptions(precision = 5, suppress = True)


#	set directory from which we work
cur_dir = os.getcwd()



#	Temp: set filename of our data
data_gene_expressions = cur_dir + "/leukemia.txt"
data_pathways = cur_dir + "/pathways.txt"



#	Numpy define our gene list, phenotypes and expression values
gene_list = list()
phenotype = np.array([], dtype = 'a3')



#	Define our expression values
expression_values = np.genfromtxt(data_gene_expressions, dtype='i4', delimiter = "\t")



#	Get our gene_list and phenotypes
def gene_list_table_open(file_name1, file_name2, gene_list, phenotype):
	#	Define phenotypes
	with open(file_name1, 'r') as temp_file:
		firstLine = temp_file.readlines()[0:1]
		for line in firstLine:
			line = line.rstrip('\n')
			split_line = re.split(r'\t+', line)
			phenotype = split_line[1 : len(split_line)]

	#	Define gene list
	with open(file_name1, 'r') as temp_file:
		for line in temp_file:
			line = line.rstrip('\n')
			split_line = re.split(r'\t+', line)
			gene_list.append(split_line[0])
	
	with open(file_name2, 'r') as temp_file, open(file_name2+'_clean.txt', 'w') as output_file:		
		#	Cleanup our pathways
		print('\nFollowing gene sets will be removed: \n')
		for line in temp_file:
			line = line.rstrip('\n')
			split_line = re.split(r'\t+', line)
			obsolete = "OBSOLETE"
			if ( (obsolete in split_line[1]) == True ) or ( len(split_line) - 2 < 15 ) :
				print('Removed:', split_line[0], '\t', 'value:', len(split_line) - 2)
			else:
				output_file.write(line + '\n')

	return gene_list, phenotype



#	Clean up our expression values
def gene_expression_values_cleanup(expression_values_array):
	expression_values_array = np.delete(expression_values_array, (0), axis = 0)
	expression_values_array = np.delete(expression_values_array, (0), axis = 1)
	return expression_values_array



#	Get signal to noise ratio and pair it with our genes
def get_signal_to_noise_sort(gene_list, phenotype, expression_values):

	array_out = np.array([])
	#	Function to calculate our signal to noise ratio
	def signal_to_noise(list1, list2):
		mju1 = np.mean(list1)
		mju2 = np.mean(list2)
		sigma1 = np.std(list1)
		sigma2 = np.std(list2)
		return ( (mju1 - mju2) / (sigma1 + sigma2) )

	#	Define our phenotypes position
	arr1 = []
	arr2 = []
	i = 0
	for pheno in phenotype:
		if pheno == 'ALL':
			arr1.append(i)
		else:
			arr2.append(i)
		i += 1
	arr1 = np.array(arr1)
	arr2 = np.array(arr2)

	i = 0

	#	Create gene/StNR array
	for gene in gene_list:
		StNR = signal_to_noise(expression_values[i,arr1], expression_values[i,arr2])
		#	Create a structured array of types 'gene' and 'stnr'
		#gene = np.core.defchararray.decode(gene)
		array_temp = np.array([(gene, StNR)], dtype = [('gene', '|S10'), ('stnr', 'f4')])
		if i == 0:
			array_out = array_temp
		else:
			array_out = np.append(array_out, array_temp)
		i = i+1

	#	Sort our gene/StNR array in descending manner
	array_out = np.sort(array_out, order = 'stnr')[::-1]
	return array_out



gene_list, phenotype = gene_list_table_open(data_gene_expressions, data_pathways, gene_list, phenotype)
gene_list = np.delete(gene_list, (0))
expression_values = gene_expression_values_cleanup(expression_values)
#	Clean up our gene list

StNR_array = get_signal_to_noise_sort(gene_list, phenotype, expression_values)


#---#
#gene_list = array of our genes
#phenotype = array of our phenotypes
#expression_values = 2D array of our gene expression values
#StNR_array = 2D array [gene\StNR]
#---#

def	scramble_phenotype(phenotype):
	np.random.shuffle(phenotype)

i = 0

#print(StNR_array['gene'])
#np.chararray.replace(StNR_array['gene'], np.core.defchararray.decode(StNR_array['gene']))

temp = list()
for x in np.nditer(StNR_array['gene'], order = 'C'):
	temp.append(np.core.defchararray.decode(x))

#------------------------------------------------------------------------------#

def create_ES_values(file_name_pathways):


	print('\nCalculating ES\n')
	



	#	Define our pathways and genes in an array
	gene_set_dict = dict()

	with open(file_name_pathways+'_clean.txt', 'r') as gene_set:
		for line in gene_set:
			line = line.rstrip('\n')
			split_line = re.split(r'\t+', line)
			pathway = split_line[0]
			gene_set = split_line[2 : len(split_line)]
			gene_set_dict[pathway] = gene_set

#	1000 permutations in order to create ESnull	
#	for i in range(0, 100):	

	#	Define our pathways and corresponding ES in an array

	#	Define our N	
	N = len(gene_list)

	ES_array = np.array([])
	ESnull_array = np.array([])
	gene_i = 0

	for pathway in gene_set_dict:

		#	Parse through whole data in order to calculate ES for
		#	each gene set

		#	Get ES for each gene set
		
		def get_ES(pathway, ES_array):

			ES = 0
			ES_list = []
			ES_temp = []

			#	Calculate Nr and define Nh
			#	Iterating through row by row				
			Nh = 0
			Nr = 0

			for x in range(0, len(temp)):
				gene = temp[x]
				r = abs(StNR_array['stnr'][x])
				if gene in gene_set_dict[pathway]:
					Nr += r
					Nh += 1
#			for StNR_pair in np.nditer(StNR_array, order = 'C'):
#				gene = np.core.defchararray.decode(StNR_pair['gene'])
#				r = abs(StNR_pair['stnr'])
#				if gene in gene_set_dict[pathway]:
#					Nr = Nr + r
#					Nh = Nh + 1

			#Define Pmiss
			Pmiss = 1 / (N - Nh)

			for x in range (0, len(temp)):
				gene = temp[x]
				r = abs(StNR_array['stnr'][x])
				if gene in gene_set_dict[pathway]:
					ES += (r / Nr)
				else:
					ES -= Pmiss
				ES_list.append(ES) 

#			#	Get all ES values
#			for StNR_pair in np.nditer(StNR_array, order = 'C'):
#				gene = np.core.defchararray.decode(StNR_pair['gene'])
#				r = abs(StNR_pair['stnr'])
#				if gene in gene_set_dict[pathway]:
#					ES = ES + (r / Nr)
#				else:
#					ES = ES - Pmiss
#				ES_list.append(ES)

			if abs(max(ES_list)) > abs(min(ES_list)):
				ES = max(ES_list)
			else:
				ES = min(ES_list)

			return(ES)

		ES = get_ES(pathway, ES_array)
		gene_i += 1


		
#		print(gene_set_series)
#		print(ES_series)

t = time.process_time()
create_ES_values(data_pathways)


elapsed_time = time.process_time() - t
print(elapsed_time)