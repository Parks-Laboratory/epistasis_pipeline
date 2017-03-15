# -*- coding: utf-8 -*-
'''
Retrieves Mouse Diversity Array genotypes from database in PLINK format

Usage: make_plink_inputs.py -h for list of options
'''

import pyodbc
import os
import sys
import time
import argparse

# Warn if file already exists
def warn_if_overwrite(output_fn):
	if os.path.isfile(output_fn):
		print('\tThe file \'' + output_fn + '\' already exists, and will be overwritten in 3 seconds (press Ctrl + C to prevent overwrite)')
		time.sleep(3)

def get_genotypes(strains, iids, output_fn, output_dir, server, db, table, idCol, chrCol, posCol):
	# Uses 0 as genetic distance (centimorgans) for all snps
	query_template = 'select ' + chrCol +','+ idCol +', 0 AS centimorgans,'+ posCol + ', %s ' +\
	'from ' + db +'.'+ table +\
	'order by ' + chrCol +','+ posCol

	# warn_if_overwrite(output_fn)

	# Create file for TPED
	output_path = os.path.join(output_dir, output_fn)
	outfile = open(output_path, 'w')
	c = pyodbc.connect(SERVER=server,DATABASE=db,DRIVER='{SQL Server Native Client 11.0}',Trusted_Connection='Yes')
	q = query_template % ', '.join(['[%s]' % x for x in strains])

	t0 = time.clock()	# see how long query took
	res = c.execute(q)
	print('\tQuery completed in %.3f minutes' % ((time.clock()-t0)/60) )

	tfam = 0
	linebuffer = []
	# maybe retrieve all rows, then print at end
	for row in res:
		# generate .tfam file for PLINK
		if tfam == 0:
			# sanitize strain names
			colnames = [t[0] for t in row.cursor_description]
			colnames = [x.replace('/', '.').replace(' ', '.') for x in colnames]

			tfam_output_path = output_path.replace('.tped', '.tfam')
			# warn_if_overwrite(tfam_output_path)
			tfam_outfile = open(tfam_output_path, 'w')
			# accesses list of strains by referring to table's column names
			for i, fid in enumerate(colnames[4:]):
				if iids is None:
					iid = (i+1)
				else:
					iid = iids[i].replace('/', '.').replace(' ', '.')
				# sets most fields in file to "missing"
				tfam_outfile.write('\t'.join(map(str, [fid, iid, 0, 0, 0, -9]))+'\n' )
			tfam_outfile.close()
			tfam = 1

		linebuffer.append('\t'.join(map(str, row)))
		# print 50000 lines at a time
		if len(linebuffer) == 50000:
			outfile.write('\n'.join(linebuffer)+'\n' )
			outfile.flush()
			linebuffer = []

	# flush lines if any are left over
	if linebuffer:
		outfile.write('\n'.join(linebuffer) )

	c.close()
	outfile.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-out', required=False, default='out_'+''.join(str(x) for x in time.gmtime())+'/plink_input',
		help='name of file-prefix to use when storing results. Creates .')
	parser.add_argument('-strains', required=True, action='append',
		help='name of file w/ column for strains (ids column optional)')
	parser.add_argument('-server', required=False, default='PARKSLAB',
		help='name of SQL server containing genotype database')
	parser.add_argument('-db', required=True,
		help='name of SQL database containing genotype tables/views')
	parser.add_argument('-table', required=True,
		help='name of SQL table/view containing genotypes for strains \
		in PLINK format')
	parser.add_argument('-idCol', required=False, default='rsID',
		help='name of column containing marker identifiers')
	parser.add_argument('-chrCol', required=False, default='snp_chr',
		help='name of column containing marker chromosome labels')
	parser.add_argument('-posCol', required=False, default='snp_bp_mm10',
		help='name of column containing marker genetic distance')
	args = parser.parse_args()
	for filename in args.strains:
		# print ('\tBuilding PLINK inputs from', filename)
		f = open(filename)
		# strip() removes whitespace from beginning/end of linebuffer
		# split() returns list of words in string, parsed using parameter char.
		lines = [x.strip().split('\t') for x in f.readlines() if x.strip()]
		f.close()

		try:
			strains, iids = zip(*lines)
		except ValueError:
			# print('\tError:',os.path.basename(sys.argv[1]), 'must have format <strain> TAB <IID>')
			# print('\tThe IIDS will be auto-generated.')
			iids = None
			# Convert lines to strings
			strains =[]
			for x in lines:
				strains.append(''.join(x))

		output_dir = os.path.split(args.out)[0]
		if( not os.path.isdir(output_dir) ):
			os.mkdir(output_dir)
		output_fn = os.path.split(args.out)[1]+'.tped'


		get_genotypes(strains, iids, output_fn, output_dir, args.server, args.db, args.table, args.idCol, args.chrCol, args.posCol)
