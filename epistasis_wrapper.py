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
fastlmmc = './fastlmmc'
plink = './plink'

species_chroms = {'human':24, 'mouse':21}
# ignore Y chromosome
species_chroms = {'human':23, 'mouse':20}

# filter out SNPs with MAF < 5%, missing genotype frequency > 10%; convert to binary format
subset_tped_file = "sed -n '%(snp_range1)sp' %(dataset)s.tped > sub%(dataset)s.tped"
subset_tped_file2 = "sed -n '%(snp_range2)sp' %(dataset)s.tped >> sub%(dataset)s.tped"
subset_tped_file3 = "cp %(dataset)s.tfam sub%(dataset)s.tfam"
make_bed_cmd = '%(plink_location)s --tfile sub%(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out %(dataset)s %(plink_species)s'
# make_bed_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --snps %(snp_range)s --make-bed --out %(dataset)s %(plink_species)s'
make_bed_all_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out all_%(dataset)s %(plink_species)s'

# calculate MAF and missing genotype frequency
# make_maf_cmd = '%(plink_location)s --bfile %(dataset)s --allow-no-sex --snps %(snp_range)s --out %(dataset)s %(plink_species)s --freq'
# make_missing_cmd = '%(plink_location)s --bfile %(dataset)s --allow-no-sex --snps %(snp_range)s --out %(dataset)s %(plink_species)s  --missing'
# clean_cmd = 'for f in %(dataset)s.{nosex,log}; do if [ -e $f ]; then rm $f; fi; done'


# function to make sure that fids and iids match across files
def check_fids_iids(prefix):
    def get_fids_iids(fn, skip=0):
        f = open(fn)
        for i in range(skip):
            f.readline()
        #fid_iid = [tuple(x.strip().split('\t')[:2]) for x in f.readlines()]
        fid_iid = [tuple(x.strip().split()[:2]) for x in f.readlines()]
        f.close()
        return fid_iid

    fid_iid_tfam, fid_iid_pheno = map(get_fids_iids, ['%s.tfam' % prefix, '%s.pheno.txt' % prefix], [0, 1])
    if fid_iid_tfam == fid_iid_pheno:
        return True
    elif set(fid_iid_tfam) == set(fid_iid_pheno):
        # FID/IID pairs can be in different order, but warn in case this is not intended
        print 'FID/IID pairs in %s.pheno.txt are not in the same order as in %s.tfam' % (prefix, prefix)
        return True
    else:
        missing = set(fid_iid_tfam).difference(fid_iid_pheno)
        extra = set(fid_iid_pheno).difference(fid_iid_tfam)

        if missing:
            # OK to have more individuals in .tfam file than .pheno.txt file
            print 'FID/IID pairs in %s.tfam but not in %s.pheno.txt:' % (prefix, prefix)
            print '\n'.join(['%s\t%s' % x for x in sorted(missing)])
            return True

        if extra:
            print 'FID/IID pairs in %s.pheno.txt but not in %s.tfam:' % (prefix, prefix)
            print '\n'.join(['%s\t%s' % x for x in sorted(extra)])

    return False


# function to run plink
def populate_available(dataset, species, snp_index):
  
    tped, tfam, pheno, covar = ['%s%s' % (dataset, suffix) for suffix in ['.tped', '.tfam', '.pheno.txt', '.covar.txt']]
    plink_species = ['', '--%s' % species][species != 'human']

    # run plink commands
    plink_location = plink
    snp_range = open('snp_combos_%s.txt' % dataset).readlines()[snp_index].strip().split(';')
    snp_range1 = snp_range[0]
    snp_range2 = snp_range[1]
    for cmd_template in (subset_tped_file, subset_tped_file2, subset_tped_file3, make_bed_cmd, make_bed_all_cmd):
        cmd = cmd_template % locals()
        print cmd
        subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)
        
    chroms = map(str, range(1, species_chroms[species] + 1))
    n_snps = int( subprocess.Popen(['wc', '-l', '%s.bim' % dataset], stdout=subprocess.PIPE).communicate()[0].split()[0] )
    sys.stdout.flush()

    if dataset in pheno:
        f = open('%s.pheno.txt' % dataset)
        # skip FID and IID columns
        headers = f.readline().strip().split('\t')[2:]
        pheno_data = f.readlines()
        f.close()

        # if FID/IID pairs match between .tfam and .pheno.txt, proceed
        # if not, go to next file
        if check_fids_iids(dataset):
            # do nothing
            print ''
        else:
            # check that phenotype names don't have spaces
            problematic = [x for x in headers if ' ' in x]
    
            # check that phenotypes exist and have unique names
            if not problematic and headers and len(set(headers)) == len(headers):
                # create phenotype key file for later reference
                f = open('pheno.key.txt', 'w')
                fmt = '%%0%sd' % int(ceil(log10(len(headers))))
                f.write( '\n'.join(['%s\t%s' % (fmt % i, phen) for i, phen in enumerate(headers)]) )
                f.close()
    
                params = {'n_pheno': len(headers),
                          'n_indivs': len(pheno_data),
                          'n_snps': n_snps}
                                                    
                if dataset in covar:
                    params['covar'] = '%s.covar.txt' % dataset
                else:
                    params['covar'] = None
    
            elif len(set(headers)) < len(headers):
                duplicates = sorted([x for x in set(headers) if headers.count(x) > 1])
                print 'Duplicated phenotype names found in %s.pheno.txt:' % dataset
                print '\n'.join(map(str, duplicates))
            elif problematic:
                print 'Spaces exist in the following phenotype names:'
                print '\n'.join(sorted(problematic))
                print "cat <(head -1 /%(dataset)s.pheno.txt | sed 's/ //g') <(tail -n+2 /%(dataset)s.pheno.txt) > tmp.txt && mv -f tmp.txt %(dataset)s.pheno.txt" % locals()
            else:
                print 'No phenotypes found in %s.pheno.txt!' % dataset


# run fastlmmc
def run_fastlmmc(dataset, output_dir, snp_index, covFile=None, species='mouse', maxthreads=1, featsel=False, exclude=False, condition=None):
    
    # commands from fastlmmc:
    # maxthreads
    # condition
    # exclude by position
    
    # if condition:
    #     condition = '-SnpId1 %s' % condition[0]
    # else:
    #     condition = ''

    bfile = dataset
    snp_reader = Bed(bfile)
    all_snp_reader = Bed('all_%s' % bfile)
    pheno = '%s.pheno.txt' % dataset

    v = globals()
    chroms = map(str, range(1, species_chroms[species] + 1))
    v.update(locals())
    

    # epistasis on all snps
    if covFile:
          
        epistasis(snp_reader, pheno, G0=all_snp_reader, covar=covFile, runner=LocalInParts(0,3,mkl_num_threads=24))
        epistasis(snp_reader, pheno, G0=all_snp_reader, covar=covFile, runner=LocalInParts(1,3,mkl_num_threads=24))
        epistasis(snp_reader, pheno, G0=all_snp_reader, covar=covFile, runner=LocalInParts(2,3,mkl_num_threads=24))
        df = epistasis(snp_reader, pheno, G0=all_snp_reader, covar=covFile, runner=LocalInParts(3,3,mkl_num_threads=24))

    else:
          
        epistasis(snp_reader, pheno, G0=all_snp_reader, runner=LocalInParts(0,3,mkl_num_threads=24))
        epistasis(snp_reader, pheno, G0=all_snp_reader, runner=LocalInParts(1,3,mkl_num_threads=24))
        epistasis(snp_reader, pheno, G0=all_snp_reader, runner=LocalInParts(2,3,mkl_num_threads=24))
        df = epistasis(snp_reader, pheno, G0=all_snp_reader, runner=LocalInParts(3,3,mkl_num_threads=24))


    # format outputs
    final = df.loc[:, ['SNP0', 'Chr0', 'ChrPos0', 'SNP1', 'Chr1', 'ChrPos1', 'PValue']]
    final.columns = ['SNP1', 'CHR1', 'BP1', 'SNP2', 'CHR2', 'BP2', 'P']
    final = final[final['P'] <= 0.00001] 

    # output to csv 
    v.update(locals())
    final.to_csv('%(output_dir)s/%(dataset)s_%(snp_index)s.gwas' % v, sep='\t', index=False)
    


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    #parser.set_usage('''%prog [options] dataset''')
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
    parser.add_argument('dataset', help='dataset to run', action='store')
    parser.add_argument('snp_index', help= 'phenotype index', action='store')

    args = parser.parse_args()

    species = args.species
    covFile = args.covFile
    featsel = args.featsel
    exclude = args.exclude
    maxthreads = args.maxthreads
    condition = args.condition
    debug = args.debug
    dataset = args.dataset
    output_dir = root
    snp_index = int( args.snp_index ) 

    if debug:
        print('args:')
        pprint(args)

    populate_available(dataset, species, snp_index)
    run_fastlmmc(dataset, output_dir, snp_index, covFile, species, maxthreads, featsel=featsel, exclude=exclude, condition=condition)


