#!/usr/bin/env python

import sys
import os
import subprocess
import csv
import struct
from math import ceil, log10
from pprint import pprint
import logging
import textwrap
from time import time
import numpy as np
from fastlmm.association import epistasis
from fastlmm.association import single_snp
from fastlmm.util.runner import LocalInParts
import pysnptools.util
from pysnptools.util.pheno import loadOnePhen
from pysnptools.snpreader import Bed

root = os.path.split(os.path.realpath(sys.argv[0]))[0]

species_chroms = {'human':24, 'mouse':21}
# ignore Y chromosome
species_chroms = {'human':23, 'mouse':20}
p_value_threshold = 0.0005

final_columns = ['SNP0','SNP1','PValue']

# run fastlmmc
def run_fastlmmc(dataset, output_dir, process_id, group_size, covFile=None, species='mouse', maxthreads=1, featsel=False, exclude=False, condition=None):

	# commands from fastlmmc:
	# maxthreads
	# condition
	# exclude by position

	# if condition:
	#	 condition = '-SnpId1 %s' % condition[0]
	# else:
	#	 condition = '

	bfile = dataset
	filtered_snp_reader = Bed('%s.FILTERED' % bfile)
	full_snp_reader = Bed('%s.FULL' % bfile)
	pheno = '%s.pheno.txt' % dataset

	v = globals()
	chroms = map(str, range(1, species_chroms[species] + 1))
	v.update(locals())

	n = len(filtered_snp_reader.sid)

	# case checking
	if(group_size > n):
		print("trying to group more than the number of existing snps:\nprogram ended!")
		exit(1)
	if(group_size == 0):
		print("grouping size is 0:\nprogram ended!")
		exit(2)
	groupNum = (n//group_size)
	if (n % group_size !=0):
		groupNum += 1
	#print("group_num: " + str(groupNum))
	if(groupNum < 2):
		print("group number should be at least two, please decrease the size of snps in each group")
		exit(3)
	th = groupNum - 1
	rest = process_id + 1
	base = 0
	hetero_num = groupNum*(groupNum - 1)//2
	maxjobnum = hetero_num
	homo_num = (groupNum //2) if (groupNum % 2 == 0) else (groupNum//2 + 1)
	maxjobnum += homo_num
	if(process_id >= maxjobnum):
		print("job number exceeds the total number of jobs that epstasis could do")
		exit(1)

	list_1_idx_start = 0
	list_1_idx_end = 0
	list_2_idx_start = 0
	list_2_idx_end = 0
	single_homo = False;
	#print("hetero_num: %s, homo_num: %s, maxjobnum: %s" %(hetero_num, homo_num,maxjobnum))
	if(process_id < hetero_num):
		while(rest > th):
			rest -= th
			th -= 1
			base += 1

		list_1_idx_start = group_size*base
		list_1_idx_end = list_1_idx_start + group_size

		list_2_idx_start = group_size*(base + rest)
		if((base + rest) == (groupNum -1)):
			 list_2_idx_end = n - 1;
		else:
			list_2_idx_end = list_2_idx_start + group_size
		#print('(' + str(base) + ',' + str(base + rest) + ')')

	# homogenous computing: same group
	else:
		offset = process_id - hetero_num

		offset *= 2
		list_1_idx_start = group_size * offset

		if(offset == groupNum - 1):
			# last homo with only one group
			#list_1_idx_start =
			list_1_idx_end = n;
			single_homo = True

		else:
			list_1_idx_end = list_1_idx_start + group_size
			list_2_idx_start = list_1_idx_end
			list_2_idx_end = min(list_2_idx_start + group_size, n)

		#print('(' + str(offset) + ',' + str(offset) + ')' + '(' + str(offset + 1) + ',' +\
		#str(offset + 1) + ')')


	# epistasis on all snps
	df = None
	df2 = None
	if covFile:
		if (process_id < hetero_num):
			df = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, covar=covFile, sid_list_0=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end], sid_list_1=filtered_snp_reader.sid[list_2_idx_start:list_2_idx_end])
		else:
			if(single_homo):
				df = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, covar=covFile, sid_list_0=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end], sid_list_1=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end])
			else:
				df = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, covar=covFile, sid_list_0=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end], sid_list_1=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end])
				df2 = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, covar=covFile, sid_list_0=filtered_snp_reader.sid[list_2_idx_start:list_2_idx_end], sid_list_1=filtered_snp_reader.sid[list_2_idx_start:list_2_idx_end])
	else:
		if (process_id < hetero_num):
			df = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, sid_list_0=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end], sid_list_1=filtered_snp_reader.sid[list_2_idx_start:list_2_idx_end])
		else:
			if(single_homo):
				df = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, sid_list_0=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end], sid_list_1=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end])
			else:
				df = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, sid_list_0=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end], sid_list_1=filtered_snp_reader.sid[list_1_idx_start:list_1_idx_end])
				df2 = epistasis(filtered_snp_reader, pheno, G0=full_snp_reader, sid_list_0=filtered_snp_reader.sid[list_2_idx_start:list_2_idx_end], sid_list_1=filtered_snp_reader.sid[list_2_idx_start:list_2_idx_end])

	def format_results(df, final_columns, threshold):
		final = df.loc[:, final_columns]
		final = final[final['PValue'] <= threshold]
		return final

	v.update(locals())
	# output to csv
	final = format_results(df, final_columns, p_value_threshold)
	final.to_csv('%(output_dir)s/%(dataset)s_%(process_id)s.gwas' % v, sep='\t', index=False)
	if(df2 is not None):
		final = format_results(df2, final_columns, p_value_threshold)
		final.to_csv('%(output_dir)s/%(dataset)s_%(process_id)s.gwas' % v, mode = 'a', sep='\t', index=False)

if __name__ == '__main__':
	from argparse import ArgumentParser
	parser = ArgumentParser()
	#parser.set_usage('''%prog [options] dataset''')
	parser.add_argument('dataset', help='dataset to run', action='store')
	parser.add_argument('group_size', type=int, help='number of snps in a group', action = 'store')
	parser.add_argument('process_id', help= 'phenotype index', action='store')

	parser.add_argument('-s', '--species', dest='species', help='mouse or human',
						default=None, action='store')
	parser.add_argument('-c', '--covariate_file', dest='covFile', help='use covariate file',
						default=None, action='store')
	parser.add_argument('-f', '--feature-selection', dest='featsel', help='perform feature selection',
						default=False, action='store_true')
	parser.add_argument('-e', '--excludeByPosition', dest='exclude', help='exclude SNPs within 2Mb of tested SNP from kinship matrix construction',
						default=False, action='store_true')
	parser.add_argument('--maxthreads', dest='maxthreads', help='max # of threads to use',
						default=1, choices=range(1,17), type=int, action='store')
	parser.add_argument('--condition', dest='condition', help='condition on a SNP',
						default=None, action='store', nargs=1)
	parser.add_argument('--debug', dest='debug', help='log debugging output',
						default=False, action='store_true')

	args = parser.parse_args()

	species = args.species
	covFile = args.covFile
	featsel = args.featsel
	exclude = args.exclude
	maxthreads = args.maxthreads
	condition = args.condition
	debug = args.debug
	dataset = args.dataset
	group_size = args.group_size

	output_dir = root
	process_id = int( args.process_id )

	if debug:
		print('args:')
		pprint(args)

	run_fastlmmc(dataset, output_dir, process_id, group_size, covFile, species, maxthreads, featsel=featsel, exclude=exclude, condition=condition)
