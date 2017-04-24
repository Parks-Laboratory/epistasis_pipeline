# Epistasis Pipeline

## Step 1: Filtering markers, making pheno file (epistasis_local.py)
1. `python epistasis_local.py <name of pheno file> [--maf] [--geno] [--genotype] [--plink] [--check] [--hold]`
	* **--maf, --geno specifies the maf and geno threshold**
	* **--genotype specifies generate ONLY _prefix_.tped and _prefix_.tfam file**
	* **--plink specifies ONLY run plink on existing _prefix_.tped and _prefix_.tfam files and generate _prefix_.FILTERED(FULL).bim/bed/fam**
	* **--check specifies ONLY check fids and iids match across generated files**
	* **--hold specifies not transfer required files to ./data/ folder**
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
2. `python epistasis_server.py <prefix> [options]`
	Use 'ls results/<prefix>/ | wc -l' to check the current returned file numbers
3. run the following command on local machine: 'scp -r <CONDOR_ADDRESS>:results <destination_directory_at_Parks_Lab>'

## Requirements
* CHTC account
* SQUID directory


## (Legacy instructions) Running the Epistasis Pipeline on UW-Madison Clusters

1. add files (tfam, tped, pheno, covar) to the **epistasis/data/** directory  
2. copy **fixpheno.sh** file into the **epistasis/data/** directory  
3. run the command `sh fixpheno.sh PREFIX` where PREFIX is the actual prefix for the data set  
4. move to the **epistasis/** directory  
5. run the command `python epistasis_pipeline.py PREFIX --covar` (with additional arguments as necessary)  
	Note: for a detailed list of additional arguments, run the command `python epistasis_pipeline.py -h`  
6. wait while pipeline runs  
7. transfer the directory **results/** back to the server  
8. check the **condor_out/** directory for any errors  
9. clear out the **data/** directory for next time  

