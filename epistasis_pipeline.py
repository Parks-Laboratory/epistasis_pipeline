#!/usr/bin/env python

"""
Runs Plink, Generates Submit Files and Processes Inputs/Outputs to Cluster
"""

import sys
import os
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
epistasis_script = os.path.join(prog_path, 'epistasis_wrapper.py')
merge_script = os.path.join(prog_path, 'merge_epistasis_output.py')

class Tee(object):
    def __init__(self, filename):
        self.logfile = open(filename, 'a')

    def send_output(self, s):
        sys.stderr.write('%s\n' %s)
        self.logfile.write('%s\n' %s)

    def close(self):
        self.logfile.close()

def timestamp():
    return datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%I-%S')


# write out snp combos into a file
def make_snp_combos(dataLoc, tfile_prefix):

  # extract snps from the tped file
  subprocess.call("awk {'print $2'} %(dataLoc)s/%(tfile_prefix)s.tped > %(dataLoc)s/snps.txt" % locals(), shell=True) 
  snps = [x.strip() for x in open('%(dataLoc)s/snps.txt' % locals()).readlines()]
  if snps[0] == '':
    snps.pop(0)

  # set up the looping params for generating the range of snps  
  first = range(1, len(snps), 7500)
  last = range(1, len(snps), 7500)
  last.pop(0)
  last.append(len(snps))

  # create range of snps
  x = ['%s,%s' % (first[i], last[i]-1) for i in range(len(first))]
  # x = ['%s-%s' % (snps[ first[i]], snps[last[i]-1]) for i in range(len(first))]

  # pairwise combos of snp ranges 
  all_combos = open('%(dataLoc)s/snp_combos_%(tfile_prefix)s.txt' % locals(), 'w')
  for i in x:
    for j in x:
      if i != j:
        s = '%(i)s;%(j)s' % locals()
        all_combos.write('%s\n' % s)

  all_combos.close()

def process(params, covar=False, memory=1024, tasks=None, species='mouse', maxthreads=1, featsel=False, exclude=False, condition=None):

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
    transfer_input_files = %(root)s/epistasis_%(dataset)s.sh,%(prog_path)s/,%(dataLoc)s/

    request_cpus = 24
    request_memory = %(use_memory)sMB
    request_disk = 48GB

    queue %(n_snp_combo)s
    ''').replace('\t*', '')

    exec_template = textwrap.dedent(
    '''#!/bin/bash

    # untar your Python installation
    tar -xzvf python.tar.gz

    # make sure the script will use your Python installation
    export PATH=$(pwd)/python/bin:$PATH

    # run your script
    python epistasis_wrapper.py %(covFile)s %(debug)s %(species)s %(maxthreads)s %(feature_selection)s %(exclude)s %(condition)s %(dataset)s $1 >& epistasis_wrapper.py.output.$1
    ''').replace('\t*', '')
    

    # set memory and max threads
    if memory is None:
        local_memory = 8000

    # generate output files
    condor_output = os.path.join(condor_output_root, dataset)
    if not os.path.exists(condor_output):
        os.makedirs(condor_output)

    job_output = os.path.join(job_output_root, dataset)
    if not os.path.exists(job_output):
        os.makedirs(job_output)

    if covar and params['covar'] is None:
        log.send_output('Specified --covar but no covariate file exists; ignored')

    if condition:
        condition = condition[0]

    params.update({'root': root, 
                   'dataLoc': dataLoc,
                   'dataset': dataset,
                   'job_output': job_output,
                   'condor_output': condor_output,
                   'covFile': ['', '-c %s' % params['covar']][covar and params['covar'] is not None],
                   'epistasis_script': epistasis_script,
                   'debug': ['', '--debug'][debug],
                   'prog_path':prog_path,
                   'timestamp':datetime.ctime(datetime.now()),
                   'species': '-s %s' % species,
                   'maxthreads':'--maxthreads %s' % maxthreads,
                   'feature_selection':['', '--feature-selection'][featsel],
                   'exclude':['', '--exclude'][exclude],
                   'condition': ['', '--condition %s' % condition][condition is not None],
                   'use_memory': local_memory, 
                   'n_snp_combo': len( open('%s/snp_combos_%s.txt' % (dataLoc, dataset)).readlines() )})

    maxthreads_option = ['', '-pe shared %s' % maxthreads][maxthreads > 1]

    submit_file = open( 'epistasis_%(dataset)s.sub' % params, 'w')
    submit_file.write( (submit_template % params).replace(',,', ',') )
    submit_file.close()

    exec_file = open( 'epistasis_%(dataset)s.sh' % params, 'w')
    exec_file.write( exec_template % params )
    exec_file.close()

    subprocess.call('chmod +x epistasis_%(dataset)s.sh' % params, shell = True)

    # submit jobs to condor
    condor_cluster = subprocess.Popen(['condor_submit', 'epistasis_%(dataset)s.sub' % params], stdout=subprocess.PIPE).communicate()[0]
    condor_cluster = re.search('\d{4,}', condor_cluster).group()
    print("Submitting Jobs to Cluster %s" % condor_cluster)
    log.send_output("%s was sent to cluster %s at %s" % (params['dataset'], condor_cluster, timestamp()))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='''Runs FaST-LMM on datasets found in specified location (looks in %s by default). Each dataset should have PLINK-formatted genotype (.tped, .tfam) and alternate phenotype (*.pheno.txt) files.  Optional covariate files should be named with the same prefix as the other files and end in .covar.txt .
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
    parser.add_argument('-c', '--covar', dest='covFile', help='use covariate file',
                        default=False, action='store_true')
    parser.add_argument('-s', '--species', dest='species', help='mouse or human',
                        default='mouse', action='store', choices=['human', 'mouse', 'dog', 'horse', 'cow', 'sheep'])
    parser.add_argument('-m', '--memory', dest='memory', help='amount of RAM (in megabytes) requested per job',
                        default=None, action='store', type=int)
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

    args = parser.parse_args()

    dataset = args.dataset[0]

    dataLoc = args.datadir
    job_output_root = args.outputdir
    list_data = args.list_dataset
    covFile = args.covFile
    memory = args.memory
    numeric = args.numeric 
    species = args.species.lower()
    maxthreads = args.maxthreads
    featsel = args.featsel
    exclude = args.exclude
    debug = args.debug
    tasks = args.tasks 
    condition = args.condition

    if debug:
        log = Tee('epistasis_pipeline-%s.log' % timestamp())
    else:
        log = Tee('/dev/null')

    if tasks and len(datasets) > 1:
        log.send_output('More than one dataset specified along with --tasks option; quitting')
        log.close()
        sys.exit(0)

    log.send_output('Searching for raw data in %s' % dataLoc)

    # initiate params
    params = {}
    params['covar'] = '%s.covar.txt' % dataset
    params['covar'] = None

    # make snp combos for epistasis
    make_snp_combos(dataLoc, dataset)
    
    # run on cluster
    process(params, covFile, memory, tasks, species=species, featsel=featsel, exclude=exclude, condition=condition)

    log.close()


