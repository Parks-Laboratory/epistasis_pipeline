'''
Arguments:    (use argparse)
    maf, default = 0.05
        pass to PLINK via --maf
    missing, default = 0.1
        pass to PLINK via --geno
    input file
        columns: Strain, Sex, trait
        columns: Strain, Sex, Covar, trait1, trait2, ... (alternative)
            if Covar only contains NA, or missing, then don't create covar file
        (use filename as prefix for all generated files)


Outputs:
    *.FILTERED.bim, *.FILTERED.bed, *.FILTERED.fam
    *.FULL.bim, *.FULL.bed, *.FULL.fam
    *.pheno.txt (altered)

Goals of script:
0) make *.tped, *.tfam
    call make_plink_inputs.py
1) make *pheno.txt
2) call fix_pheno    ... replaces missing values with -9 (by convention)
    see fix_pheno
3) check fids/iids ...compares *.tfam with *.pheno.txt, makes sure 1st column same for both
    see check_fids_iids()
4) Filter snps, make beds (populate_available)    ...
    see populate_available()

call make_plink_input.py
optional arguments:
  -h, --help        show this help message and exit
  -out OUT          name of file-prefix to use when storing results. Creates .
  -strains STRAINS  name of file w/ column for strains (ids column optional)
  -server SERVER    name of SQL server containing genotype database
  -db DB            name of SQL database containing genotype tables/views
  -table TABLE      name of SQL table/view containing genotypes for strains in
                    PLINK format
  -idCol IDCOL      name of column containing marker identifiers
  -chrCol CHRCOL    name of column containing marker chromosome labels
  -posCol POSCOL    name of column containing marker genetic distance

  make_plink_inputs.py -strains %OUTPUT_DIR%/strains_list.txt -db HMDP -table [dbo].[genotype_calls_plink_format] -out %OUTPUT_DIR%/plink_input -idCol rsID -chrCol snp_chr -posCol snp_bp_mm10
'''

import argparse
import subprocess
import sys
from make_plink_inputs import get_genotypes
from math import ceil, log10

parser = argparse.ArgumentParser()
parser.add_argument('--maf', action='store', default = 0.05)
parser.add_argument('--geno', action='store', default = 0.1)
args = parser.parse_args()
maf = args.maf
geno = args.geno

make_bed_cmd = '%(plink_location)s --tfile sub%(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out %(dataset)s %(plink_species)s'
# make_bed_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --snps %(snp_range)s --make-bed --out %(dataset)s %(plink_species)s'
make_bed_all_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out all_%(dataset)s %(plink_species)s'

input_file = "EPISTASIS_TEST_TRAIT.TXT"
strains =  ([x.split('\t')[1] for x in open(input_file).readlines()][1:])
pheno_strains = [strain.replace('/', '.').replace(' ', '.') for strain in strains]
traits = ([x.split('\t')[-1] for x in open(input_file).readlines()][1:])

f = open(input_file.split(".")[0] + ".pheno.txt", "w")
for i in range(0, len(strains)):
    f.write( pheno_strains[i] + "\t")
    f.write( str(i) + "\t")
    f.write( traits[i])
f.close()
output_file_name = input_file.split(".")[0] +".tped"
#get_genotypes(strains , output_fn= output_file_name, output_dir= "", db= "HMDP",
 #             view="[dbo].[genotype_calls_plink_format]")
# Change this later
subset_tped_file = ""
subset_tped_file2 = ""
subset_tped_file3 = ""
species_chroms = []


# function to run plink (used to be in epistasis_wrapper.py)
def populate_available(dataset, species, snp_index,maf,geno):

    tped, tfam, pheno, covar = ['%s%s' % (dataset, suffix) for suffix in ['.tped', '.tfam', '.pheno.txt', '.covar.txt']]
    plink_species = ['', '--%s' % species][species != 'human']

    # run plink commands
    plink_location = 'plink'
    snp_range = open('snp_combos_%s.txt' % dataset).readlines()[snp_index].strip().split(';')
    snp_range1 = snp_range[0]
    snp_range2 = snp_range[1]
    for cmd_template in (subset_tped_file, subset_tped_file2, subset_tped_file3, make_bed_cmd, make_bed_all_cmd):
        cmd = cmd_template % locals()
        print (cmd)
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
            # print ''
            print ("")
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
                print ('Duplicated phenotype names found in %s.pheno.txt:' % dataset)
                print ('\n'.join(map(str, duplicates)))
            elif problematic:
                print ('Spaces exist in the following phenotype names:')
                print ('\n'.join(sorted(problematic)))
                print ("cat <(head -1 /%(dataset)s.pheno.txt | sed 's/ //g') <(tail -n+2 /%(dataset)s.pheno.txt) > tmp.txt && mv -f tmp.txt %(dataset)s.pheno.txt" % locals())
            else:
                print ('No phenotypes found in %s.pheno.txt!' % dataset)


# function to make sure that fids and iids match across files (used to be in epistasis_wrapper.py)
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
        print ('FID/IID pairs in %s.pheno.txt are not in the same order as in %s.tfam' % (prefix, prefix))
        return True
    else:
        missing = set(fid_iid_tfam).difference(fid_iid_pheno)
        extra = set(fid_iid_pheno).difference(fid_iid_tfam)

        if missing:
            # OK to have more individuals in .tfam file than .pheno.txt file
            print ('FID/IID pairs in %s.tfam but not in %s.pheno.txt:' % (prefix, prefix))
            print ('\n'.join(['%s\t%s' % x for x in sorted(missing)]))
            return True

        if extra:
            print ('FID/IID pairs in %s.pheno.txt but not in %s.tfam:' % (prefix, prefix))
            print ('\n'.join(['%s\t%s' % x for x in sorted(extra)]))

    return False


# call this
fixpheno='''\
#!/bin/sh


while (("$#")); do
    if [ -e $1.pheno.txt ]; then
        echo Fixing $1.pheno.txt;
                head -1 $1.pheno.txt > tmp.pheno.txt;
                tail -n+2 $1.pheno.txt | sed -r 's/ /\./g;s/\//\./g;s/(NULL|NA|#NUM!|-Inf|Inf)/-9/g' >> tmp.pheno.txt;
                mv -f tmp.pheno.txt $1.pheno.txt;
        fi;
        if [ -e $1.covar.txt ]; then
                echo Fixing $1.covar.txt;
                sed -i 's/ /\./g;s/\//\./g;s/(NULL|NA|#NUM!|-Inf|Inf)/-9/g' $1.covar.txt;
        fi;
        shift;
done'''




'''
dataset is file prefix
species is 'mouse' or 'human'
snp_index was a CHTC process number
'''
dataset = ""
species = ""
snp_index  = ""
# populate_available(dataset, species, snp_index,maf, geno)
