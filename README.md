# Epistasis Pipeline

## Step 1: Filtering markers, making pheno file (epistasis_local.py)
1. Manually generate input file with one row per individual, and the following tab-separated columns: 
	
	Strain	Sex	_trait_1_	_trait_2_ ... _trait_n_ 	Covar
	
1. `python epistasis_local.py <name of input file> [--maf] [--geno] [--genotype] [--plink] [--check] [--hold]`
	* --maf, --geno specifies the maf and geno threshold
	* --genotype specifies generate ONLY _prefix_.tped and _prefix_.tfam file
	* --plink specifies ONLY run plink on existing _prefix_.tped and _prefix_.tfam files and generate _prefix_.FILTERED(FULL).bim/bed/fam**
	* --check specifies ONLY check fids and iids match across generated files
	* --hold specifies not transfer required files to ./data/ folder
## Step 2: Preparing files for transfer, submit jobs (epistasis_server.py)
1. scp **data/**, **scripts/**, **epistasis_server.py** to submit server
	* **data/** contains 
		**_prefix_.FILTERED.bim**, 
		**_prefix_.FILTERED.bed**, 
		**_prefix_.FILTERED.fam**, 
		**_prefix_.FULL.bim**, 
		**_prefix_.FULL.bed**, 
		**_prefix_.FULL.fam**, and 
		**_prefix_.pheno.txt**
	* **scripts/** contains 
		**python.tar.gz** (portable Python 2.7 installation), 
		**atlas.tar.gz** (ATLAS library), and 	
		**epistasis_node.py**
1. `python epistasis_server.py <prefix> [options]`
	Use `ls results/<prefix> | wc -l` to check the current returned file numbers
1. run the following command on local machine: `scp -r <CONDOR_ADDRESS>:results <destination_directory_at_Parks_Lab>`

## Requirements
* CHTC account
* SQUID directory (CHTC uses uses this to transfer large files)

