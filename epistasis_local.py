'''
Goals of script:
1) call fix_pheno?
2) Filter snps, make beds
3) check fids/iids
4) make tarball containing all files and directories


'''



make_bed_cmd = '%(plink_location)s --tfile sub%(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out %(dataset)s %(plink_species)s'
# make_bed_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --snps %(snp_range)s --make-bed --out %(dataset)s %(plink_species)s'
make_bed_all_cmd = '%(plink_location)s --tfile %(dataset)s --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out all_%(dataset)s %(plink_species)s'


# function to run plink (used to be in epistasis_wrapper.py)
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
			print ('')
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


# function to make sure that fids and iids match across files (used to be in epistasis_pipeline.py)
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





'''
dataset is file prefix
species is 'mouse' or 'human'
snp_index was a CHTC process number
'''
# populate_available(dataset, species, snp_index)
