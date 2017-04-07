#!/usr/bin/env python

"""
Runs Plink, Generates Submit Files and Processes Inputs/Outputs to Cluster
"""
from __future__ import division
import sys
import os
import pwd
import stat
import re
import subprocess
from math import log10, ceil
import textwrap
import operator
from datetime import datetime
import time
import re



debug = False

# input/output folders
root = os.path.split(os.path.realpath(sys.argv[0]))[0]
dataLoc = os.path.join(root, 'data')
condor_output_root = os.path.join(root, 'condor_out')
job_output_root = os.path.join(root, 'results')

# script locations
prog_path = os.path.join(root, 'scripts')
epistasis_script = os.path.join(prog_path, 'epistasis_node.py')

# filename formats
FULL_DATASET = '.FULL'
FILTERED_DATASET = '.FILTERED'

class Tee(object):
	def __init__(self, filename):
		self.logfile = open(filename, 'a')

	def send_output(self, s):
		sys.stderr.write('%s\n' %s)
		self.logfile.write('%s\n' %s)

	def close(self):
		self.logfile.close()

def timestamp():
	return datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%M-%S')

def process(params):

	os.chdir(root)

	# submit and executable files
	submit_template = textwrap.dedent(
	'''# Epistasis Submit File

	universe = vanilla
	log = %(condor_output)s/epistasis_$(Cluster).log
	error = %(condor_output)s/epistasis_$(Cluster)_$(Process).err

	InitialDir = %(root)s/results/%(dataset)s
	executable = %(root)s/epistasis_%(dataset)s.sh
	arguments = $(Process)
	output = %(condor_output)s/epistasis_$(Cluster)_$(Process).out

	should_transfer_files = YES
	when_to_transfer_output = ON_EXIT
	transfer_input_files = http://proxy.chtc.wisc.edu/SQUID/%(username)s/%(squid_zip)s

	request_cpus = 1
	request_memory = %(use_memory)sMB
	request_disk = 2GB

	# set the interval for releasing the job if failed
	periodic_release = (CurrentTime - EnteredCurrentStatus > 600)

	# requirements = (Target.PoolName =!= "CHTC")
	+wantGlidein = true
	+wantFlocking = true

	queue %(num_jobs)s
	''').replace('\t*', '')


	exec_template = textwrap.dedent(
	'''#!/bin/bash

	# if script fails before getting to python, make sure this file exists
	echo > epistasis_node.py.output.$1

	# untar your files sent along by SQUID
	tar -xzvf %(squid_zip)s

	# untar your Python installation
	tar -xzvf %(python_installation)s

	# untar ATLAS linear algebra library
	tar -xzvf %(atlas_installation)s

	# make sure the script will use your Python installation
	export PATH=$(pwd)/python/bin:$PATH

	# make sure script can find ATLAS library
	export LD_LIBRARY_PATH=$(pwd)/atlas

	# run your script
	python epistasis_node.py %(dataset)s %(group_size)s $1 %(covFile)s %(debug)s %(species)s %(maxthreads)s %(feature_selection)s %(exclude)s %(condition)s &>> epistasis_node.py.output.$1

	# Keep job output only if job FAILS (for debugging/so it can be re-run),
	# otherwise, delete the output file
	if [ $? == 0 ]; then
		rm epistasis_node.py.output.$1
	fi

	rm -r -f *.bed *.bim *.fam *.py *.pyc *.tar.gz *.txt python
	''').replace('\t*', '')



	# generate output files
	if not os.path.exists(params['condor_output']):
		os.makedirs(params['condor_output'])

	if not os.path.exists(params['job_output']):
		os.makedirs(params['job_output'])

	submit_file = open( 'epistasis_%(dataset)s.sub' % params, 'w')
	submit_file.write( (submit_template % params).replace(',,', ',') )
	submit_file.close()

	exec_file = open( 'epistasis_%(dataset)s.sh' % params, 'w')
	exec_file.write( exec_template % params )
	exec_file.close()

	# give script permission to execute
	subprocess.call('chmod +x epistasis_%(dataset)s.sh' % params, shell = True)

	# add large files to a tar archive file, which will be sent over by SQUID
	subprocess.call('tar -cf %(squid_archive)s -C %(dataLoc)s/ .' % params, shell = True)
	subprocess.call('tar -f %(squid_archive)s -C %(prog_path)s --append .' % params, shell = True)
	# compress archive file
	subprocess.call('gzip < %(squid_archive)s > %(squid_zip)s' % params, shell = True)
	# place compressed archive file in the user's SQUID directory
	if(subprocess.call('cp %(squid_zip)s /squid/%(username)s' % params, shell = True)):
		sys.exit('Failed to create %(squid_zip)s and copy it to squid directory' % params)

	submit_jobs(params)

def write_submission_files(params):
	pass

def write_submission_file(params, offset):
	pass

def write_shell_script(params):
	pass

def make_output_dirs(params):
	pass

def package_SQUID_files(params):
	pass

def submit_jobs(params):
	# submit jobs to condor
	condor_cluster = subprocess.Popen(['condor_submit', 'epistasis_%(dataset)s.sub' % params], stdout=subprocess.PIPE).communicate()[0]
	condor_cluster = re.search('\d{4,}', condor_cluster).group()
	print("Submitting Jobs to Cluster %s" % condor_cluster)
	log.send_output("%s was sent to cluster %s at %s" % (params['dataset'], condor_cluster, timestamp()))

def num_jobs(group_size):
	num_snps = 0
	with open(os.path.join(dataLoc, dataset+FILTERED_DATASET+'.bim')) as file:
		for line in file.readlines():
			line = line.split()
			if len(line) > 1 and 'rs' in line[1]:
				num_snps += 1

	num_groups = ceil(num_snps/group_size)

	# for all groups A, B: num jobs comparing A to B, but ignoring the redundant
	# jobs comparing B to A
	num_AB_jobs = num_groups * (num_groups - 1) / 2

	# for each comparison A to A, B to B that have not been done,
	# 2 such comparisons are done per job
	# (if num_groups is odd, then one job will only compare A to A, and not B to B)
	num_AA_jobs = ceil(num_groups/2)

	return int(num_AB_jobs + num_AA_jobs)

def check_prefixes(dataloc, dataset):
	'''
	Looks in directory specified by dataloc, and ensures there are 2 sets of
	binary files. Specifically, the following must exist:
	*.FULL.bed, *.FULL.bim, *.FULL.fam
	*.FILTERED.bed, *.FILTERED.bim, *.FILTERED.fam
	where * is the same for all 6 files
	'''
	bin_files = [x for x in os.listdir(dataloc) if os.path.splitext(x)[1] in ['.bed', '.bim', '.fam']]
	full_prefixes = [x for x in bin_files if os.path.splitext(x)[0] == dataset+FULL_DATASET]
	filtered_prefixes = [x for x in bin_files if os.path.splitext(x)[0] == dataset+FILTERED_DATASET]

	# check that prefixes match up and that each extension is included exactly twice
	if not len(full_prefixes) == 3 or not len(filtered_prefixes) == 3:
		sys.exit('ERROR: two sets of .bed/.bim/.fam files could not be found in "{}" with the prefix "{}"'.format(dataloc, dataset))



if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='''Runs FaST-LMM on datasets found in specified location (looks in %s by default). Each dataset should have PLINK-formatted genotype (*.tped, *.tfam) and alternate phenotype (*.pheno.txt) files.  Optional covariate files should be named with the same prefix as the other files and end in .covar.txt .
''' % dataLoc)
	#parser.set_usage('''%(prog)s [options] [dataset1] [dataset2] ... (runs all datasets if unspecified)
	#PLINK-formatted genotype (*.tped, *.tfam) and alternate phenotype (*.pheno.txt) files should be placed in ''' + dataLoc)
	parser.add_argument('dataset', metavar='dataset', nargs=1, type=str, help='dataset(s) to process')

	parser.add_argument('-l', '--list', dest='list_dataset', help='lists datasets to process, does not do processing',
						default=False, action='store_true')
	parser.add_argument('-d', '--datadir', dest='datadir', help='specifies folder to search for raw data',
						default=dataLoc, action='store')
	parser.add_argument('-o', '--outputdir', dest='outputdir', help='specifies output folder',
						default=job_output_root, action='store')
	parser.add_argument('-c', '--covar', dest='covar', help='use covariate file',
						default=False, action='store_true')
	parser.add_argument('-s', '--species', dest='species', help='mouse or human',
						default='mouse', action='store', choices=['human', 'mouse', 'dog', 'horse', 'cow', 'sheep'])
	parser.add_argument('-m', '--memory', dest='memory', help='amount of RAM (in megabytes) requested per job',
						default=2000, action='store', type=int)
	parser.add_argument('--maxthreads', dest='maxthreads', help='maximum # of threads to use',
						default=1, action='store', choices=range(1, 17), type=int)
	parser.add_argument('-f', '--feature-selection', dest='featsel', help='perform feature selection',
						default=False, action='store_true')
	parser.add_argument('-e', '--excludeByPosition', dest='exclude', help='exclude SNPs within 2Mb of tested SNP from kinship matrix construction',
						default=False, action='store_true')
	parser.add_argument('-n', '--numeric_phenotype_id', dest='numeric', help='convert phenotype names to numbers (for safety)',
						nargs='?', default=0, const=1, type=int, action='store', choices=[0, 1, 2])
	parser.add_argument('-q', '--quiet', dest='debug', help="suppress debugging output",
						default=True, action='store_false')
	parser.add_argument('--tasks', dest='tasks', metavar='TASK', nargs='+', help='run only specified sub-tasks (specify only one dataset when using this option)', type=int)
	parser.add_argument('--condition', dest='condition', help='condition on SNP {snp_id}',
						action='store', nargs=1)
	parser.add_argument('-g', '--group_size', type=int, help='number of snps in a group', action = 'store', default=1200)

	args = parser.parse_args()

	dataset = args.dataset[0]

	dataLoc = args.datadir
	job_output_root = args.outputdir
	list_data = args.list_dataset
	covar = args.covar
	memory = args.memory
	numeric = args.numeric
	species = args.species.lower()
	maxthreads = args.maxthreads
	featsel = args.featsel
	exclude = args.exclude
	debug = args.debug
	tasks = args.tasks
	condition = args.condition
	group_size = args.group_size

	if debug:
		log = Tee('epistasis_pipeline-%s.log' % timestamp())
	else:
		log = Tee('/dev/null')

	if tasks and len(datasets) > 1:
		log.send_output('More than one dataset specified along with --tasks option; quitting')
		log.close()
		sys.exit(0)

	log.send_output('Searching for raw data in %s' % dataLoc)

	params = {}

	# initiate params
	covFile = '%s.covar.txt' % dataset
	params.update({ 'covFile': ''})
	if covar and os.path.isfile(covFile):
		params.update({ 'covFile': '-c %s' % covFile})
	elif covar:
		log.send_output('Specified --covar but no covariate file exists; ignored')

	if condition:
		condition = condition[0]

	# check_prefixes(dataLoc, dataset)
	squid_archive = 'epistasis.tar'
	params.update({'root': root,
				   'dataLoc': dataLoc,
				   'dataset': dataset,
				   'group_size': group_size,
				   'num_jobs': num_jobs(group_size),
				   'job_output': os.path.join(job_output_root, dataset),
				   'condor_output': os.path.join(condor_output_root, dataset),
				   'epistasis_script': epistasis_script,
				   'squid_archive': squid_archive,
				   'squid_zip': squid_archive + '.gz',
				   'username': pwd.getpwuid(os.getuid()).pw_name,
				   'python_installation': 'python.tar.gz',
				   'atlas_installation': 'atlas.tar.gz',
				   'debug': ['', '--debug'][debug],
				   'prog_path':prog_path,
				   'timestamp':datetime.ctime(datetime.now()),
				   'species': '-s %s' % species,
				   'maxthreads':'--maxthreads %s' % maxthreads,
				   'feature_selection':['', '--feature-selection'][featsel],
				   'exclude':['', '--exclude'][exclude],
				   'condition': ['', '--condition %s' % condition][condition is not None],
				   'use_memory': memory})

	# maxthreads_option = ['', '-pe shared %s' % maxthreads][maxthreads > 1]

	# run on cluster
	process(params)

	log.close()
