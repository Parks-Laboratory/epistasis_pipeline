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
    *.pheno.txt

Goals of script:
0) maybe prune out hetero strains
1) make *.tped, *.tfam
    call make_plink_inputs.py
2) make *pheno.txt
4) call fix_pheno    ... replaces missing values with -9 (by convention)
    see fix_pheno
8) check fids/iids ...compares *.tfam with *.pheno.txt, makes sure 1st column same for both
    see check_fids_iids()
16) Filter snps, make beds (populate_available)    ...
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
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument('file', action = 'store')
parser.add_argument('--maf', action='store', default = 0.05)
parser.add_argument('--geno', action='store', default = 0.1)
parser.add_argument('--covar', action='store_true', default=False)
'''
--genotype: select to run get genotype 
--plink: select to run plink command 
--check: select to run check IID/FID
--hold: hold the ouput file to the current directory
'''
parser.add_argument('--genotype', action='store_true', default=False)
parser.add_argument('--plink', action='store_true', default=False)
parser.add_argument('--check', action='store_true', default=False)
parser.add_argument('--hold', action='store_true', default=False)

args = parser.parse_args()
if ((not args.plink) and (not args.check) and (not args.genotype)):
    plink = check = genotypes = True
else:
    plink, check, genotypes = args.plink, args.check, args.genotype

hold = args.hold
maf = args.maf
geno = args.geno
covar = args.covar

make_bed_cmd = '%(plink_location)s --tfile sub%(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out %(dataset)s %(plink_species)s'
# make_bed_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --snps %(snp_range)s --make-bed --out %(dataset)s %(plink_species)s'
make_bed_all_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out all_%(dataset)s %(plink_species)s'
def convert_missing_value(str):
    grab = re.sub(pattern= " ", repl= ".", string= str)
    if (grab == "NULL" )or (grab == "NA") or (grab  == "#NUM!" )or (grab  == "-Inf") or (grab == "Inf"):
        grab = "-9"
    return grab
# print("args.file: " + args.file)
input_file = args.file
prefix = input_file.split(".")[0]
with open(input_file) as f:
    traits = f.readline().split("\t")
    last_column =traits[-1]
suffix = ".covar.txt" if (last_column.strip() == 'Covar') else ".pheno.txt"
header = ["FID", "IID", "\t".join(traits[2:])]
header = "\t".join(header)

strains =  ([x.split('\t')[0] for x in open(input_file).readlines()][1:])
pheno_strains = [strain.replace('/', '.').replace(' ', '.') for strain in strains]
traits = ([x.split('\t')[2:] for x in open(input_file).readlines()][1:])

# fixphenos and write to %s.pheno.txt %pheno_prefix
f = open(prefix + suffix, "w")
f.write(header)
for i in range(0, len(strains)):
    # replace (NULL|NA|#NUM!|-Inf|Inf) with -9
    # replace " " with "\."
    pheno_strains[i] = convert_missing_value (pheno_strains[i])
    for trait in traits[i]:
        trait = convert_missing_value (trait)
    f.write( pheno_strains[i] + "\t")
    f.write( str(i) + "\t")
    f.write( ".\t".join(traits[i]))
f.close()
print("fixed pheno and generated %s file" %(suffix))

''''
Call get_genotypes to get the .tped and .tfam file
'''
if  genotypes:
    get_genotypes(strains , output_fn= prefix, output_dir= "", db= "HMDP",
               view="[dbo].[genotype_calls_plink_format]")

print("generated  .tped and .tfam file using get_genotypes")

# function to run plink (used to be in epistasis_wrapper.py)
def populate_available(dataset, maf,geno):
    # plink_species = ['', '--%s' % species][species != 'human']

    # run plink commands
    plink_location = 'plink'
    cmd1 = "plink --tfile %s --make-bed -out %s.FULL" %(prefix, prefix)
    cmd2 = "plink --bfile %s.FULL --maf %f --geno %f --make-bed  -out %s.FILTERED" %(prefix, maf, geno, prefix)
    f = open("plink_stdout.txt", "w")
    subprocess.check_call(cmd1, shell=True, stderr=subprocess.STDOUT, stdout=f)
    subprocess.check_call(cmd2, shell=True, stderr=subprocess.STDOUT, stdout=f)
    f.close()
    #chroms = map(str, range(1, species_chroms[species] + 1))
    # n_snps = int(subprocess.Popen(['wc', '-l', '%s.bim' % dataset], stdout=subprocess.PIPE).communicate()[0].split()[0])
    sys.stdout.flush()


def check_headers(dataset):
    global suffix
    tped, tfam, pheno, covar = ['%s%s' % (dataset, all_suffix) for all_suffix in ['.tped', '.tfam', suffix, '.covar.txt']]
    if dataset in pheno:
        f = open('%s%s' % (dataset, suffix))
        # skip FID and IID columns
        headers = f.readline().strip().split('\t')[2:]
       # pheno_data = f.readlines()
        f.close()

        # if FID/IID pairs match between .tfam and suffix, proceed
        # if not, go to next file
        if check_fids_iids(dataset):
            # do nothing
            # print ''
            print ("")
        else:
            # check that phenotype names don't have spaces
            problematic = [x for x in headers if ' ' in x]
            param = {"suffix":suffix, "dataset": dataset}
            # check that phenotypes exist and have unique names
            if not problematic and headers and len(set(headers)) == len(headers):
                # create phenotype key file for later reference
                f = open('pheno.key.txt', 'w')
                fmt = '%%0%sd' % int(ceil(log10(len(headers))))
                f.write('\n'.join(['%s\t%s' % (fmt % i, phen) for i, phen in enumerate(headers)]))
                f.close()

                # params = {'n_pheno': len(headers),
                #           'n_indivs': len(pheno_data)}
                #           #\,'n_snps': n_snps}

                # if dataset in covar:
                #     params['covar'] = '%s.covar.txt' % dataset
                # else:
                #     params['covar'] = None

            elif len(set(headers)) < len(headers):
                duplicates = sorted([x for x in set(headers) if headers.count(x) > 1])
                print ('Duplicated phenotype names found in %s%s:' % (dataset, suffix))
                print ('\n'.join(map(str, duplicates)))
            elif problematic:
                print ('Spaces exist in the following phenotype names:')
                print ('\n'.join(sorted(problematic)))
                print ("cat <(head -1 /%(dataset)s%(suffix)s| sed 's/ //g') <(tail -n+2 /%(dataset)s%(suffix)s) > tmp.txt && mv -f tmp.txt %(dataset)s%(suffix)s" % param)
            else:
                print ('No phenotypes found in %(dataset)s%(suffix)s!' % param)

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
    global suffix
    fid_iid_tfam, fid_iid_pheno = map(get_fids_iids, ['%s.tfam' % prefix, '%s%s' % (prefix, suffix)], [0, 1])
    if fid_iid_tfam == fid_iid_pheno:
        return True
    elif set(fid_iid_tfam) == set(fid_iid_pheno):
        # FID/IID pairs can be in different order, but warn in case this is not intended
        print ('FID/IID pairs in %s%s are not in the same order as in %s.tfam' % ( prefix, suffix, prefix))
        return True
    else:
        missing = set(fid_iid_tfam).difference(fid_iid_pheno)
        extra = set(fid_iid_pheno).difference(fid_iid_tfam)

        if missing:
            # OK to have more individuals in .tfam file than .pheno.txt file
            print ('FID/IID pairs in %s.tfam but not in %s%s:' % (prefix, prefix, suffix))
            print ('\n'.join(['%s\t%s' % x for x in sorted(missing)]))
            return True

        if extra:
            print ('FID/IID pairs in %s%s but not in %s.tfam:' % (prefix, prefix, suffix))
            print ('\n'.join(['%s\t%s' % x for x in sorted(extra)]))

    return False

def transfer_files():
    import shutil, fnmatch
    # move data .bed, *.bim, *.fam, *.pheno.txt into new directory data
    if not os.path.exists("data"):
        os.makedirs("data")

    # filelist = [file for file in os.listdir('.') if (re.search(".*.bed|.*.bim|.*.fam|.*{0!s}".format(suffix), file))]
    filelist2 = [file for file in os.listdir('.') if (re.search("{0!s}.*\.bed|{0!s}.*\.bim|{0!s}.*\.fam|{0!s}.*{1!s}"
                                                                .format(prefix, suffix).replace("'", ""), file))]

    for file in filelist2:
        if(os.path.exists("./data/{0!s}".format(file))):
            print("file:{0!s} already exists in /data, and will be overwritten...".format(file))
            os.remove("./data/{0!s}".format(file))
        shutil.move(file, "./data", )
        print("successfully moved file:{0!s} to ./data".format(file))

    print("moved following files to directory data: \n" + "\n".join(filelist2))
if check:
    check_fids_iids(prefix)
print("finished checking the fids and iids")

if plink:
    populate_available(prefix, maf, geno)
    print ("finished generating .bed,.bim,.ped file filetering the snps specified by --maf and --geno using plink")

check_headers(prefix)
print("finsihed checking the unique values of headers")

if not hold:
    transfer_files()

'''
dataset is file prefix
species is 'mouse' or 'human'
snp_index was a CHTC process number
'''
dataset = ""
species = ""
snp_index  = ""
# populate_available(dataset, species, snp_index,maf, geno)
