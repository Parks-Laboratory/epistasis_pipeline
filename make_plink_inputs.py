# -*- coding: utf-8 -*-
'''
Retrieves Mouse Diversity Array genotypes from database in PLINK format for
the strains specified in input file. Writes .tped and .tfam files.

Get description of script parameters:		make_plink_inputs.py -h
'''

import pyodbc
import os
import sys
import time
import argparse

defaults = {
	'HMDP': {
		'table': '[dbo].[genotype_calls_plink_format]',
		'idCol': 'rsID',
		'posCol': 'snp_bp_mm10',
		'chrCol': 'snp_chr',
	},
	'DO': {
		'table': '',
		'idCol': 'snp_id',
		'posCol': 'snp_bp_mm10',
		'chrCol': 'snp_chr',
	},
}

# Warn if file already exists
def warn_if_overwrite(output_fn):
	if os.path.isfile(output_fn):
		print('\tThe file \'' + output_fn + '\' already exists, and will be overwritten in 3 seconds (press Ctrl + C to prevent overwrite)')
		time.sleep(3)

def get_genotypes(strains, output_fn, db, table=None, server=None, idCol=None, chrCol=None, posCol=None, iids=None, output_dir=None):
	'''
	Arguments:
	strains -- list of strain names
	iids -- list of IDs. If None, these will be generated
	output_fn -- prefix for generated .tped, .tfam files
	output_dir -- location where .tped, .tfam files should be placed.
	db -- SQL database
	table -- SQL table containing genotypes
	server -- SQL server
	idCol -- column in table containing marker identifiers (e.g. rsID)
	chrCol -- column in table containing marker chromosome labels (e.g. snp_chr)
	posCol -- column in table containing marker genetic distance (e.g. snp_bp_mm10)
	'''

	if table is None:
		table=defaults[db]['table']
	if server is None:
		server='PARKSLAB'
	if idCol is None:
		idCol=defaults[db]['idCol']
	if chrCol is None:
		chrCol=defaults[db]['chrCol']
	if posCol is None:
		posCol=defaults[db]['posCol']
	if output_dir is None:
		output_dir = ''
	elif output_dir and not os.path.isdir(output_dir):
		os.mkdir(output_dir)

	query_template = 'select ' + chrCol +','+ idCol +', 0 AS centimorgans,'+ posCol + ', %s ' +\
	'from ' + db +'.'+ table +\
	'order by ' + chrCol +','+ posCol

	# warn_if_overwrite(output_fn)

	output_fn += '.tped'
	# Create file for TPED
	output_path = os.path.join(output_dir, output_fn)
	outfile = open(output_path, 'w')
	c = pyodbc.connect(SERVER=server,DATABASE=db,DRIVER='{SQL Server Native Client 11.0}',Trusted_Connection='Yes')
	q = query_template % ', '.join(['[%s]' % x for x in strains])

	# t0 = time.clock()	# see how long query took
	res = c.execute(q)
	# print('\tQuery completed in %.3f minutes' % ((time.clock()-t0)/60) )

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
					# PLINK requires IID to be > 0
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
	parser = argparse.ArgumentParser(description='''Given a (tab-separated) input file
		containing a column of strains and (optionally) a column of IDs,
		query an SQL table for PLINK-formatted
		genotype data for only the strains specified, and
		use this data to write .tped and .tfam files. (See plink_table.sql for format of SQL table)''')

	parser.add_argument(dest='strains', action='append',
		help='file w/ tab-separated format: [strains column] [IDs columns (optional)]')
	parser.add_argument(dest='db',
		help='SQL database containing genotype table (e.g. HMDP)')

	parser.add_argument('-table',
		help='SQL table containing genotypes for strains \
		in PLINK format (e.g. [dbo].[genotype_calls_plink_format])')
	parser.add_argument('-output_dir', required=False,
		help='directory in which to store results.')
	parser.add_argument('-server', required=False,
		help='SQL server containing genotype database (e.g. PARKSLAB)')
	parser.add_argument('-idCol', required=False,
		help='column in table containing marker identifiers (e.g. rsID)')
	parser.add_argument('-chrCol', required=False,
		help='column in table containing marker chromosome labels (e.g. snp_chr)')
	parser.add_argument('-posCol', required=False,
		help='column in table containing marker genetic distance (e.g. snp_bp_mm10)')
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

		get_genotypes(strains=strains, output_fn=os.path.splitext(filename)[0], db=args.db, table=args.table, server=args.server, idCol=args.idCol, chrCol=args.chrCol, posCol=args.posCol, iids=iids, output_dir=args.output_dir)
