# Epistasis Pipeline
[Documentation](http://microsoftgenomics.github.io/FaST-LMM/#epistasis)

## Step 1: Filtering markers, making pheno file (epistasis_local.py)
1. Manually generate input file with one row per individual, and the following tab-separated columns:

	Strain	Sex	_trait_ 	Covar
	
	Note: the "Covar" column is optional

1. `python epistasis_local.py <name of input file> [--maf] [--geno] [--genotype] [--plink] [--check] [--hold]`
	* --maf, --geno specifies the maf and geno threshold
	* --genotype specifies generate ONLY _prefix_.tped and _prefix_.tfam file
	* --plink specifies ONLY run plink on existing _prefix_.tped and _prefix_.tfam files and generate _prefix_.FILTERED(FULL).bim/bed/fam**
	* --check specifies ONLY check fids and iids match across generated files
	* --hold specifies not transfer required files to ./data/ folder
## Step 2: Preparing files for transfer, submit jobs (epistasis_submit.py)
1. scp **data/**, **scripts/**, **epistasis_submit.py** to submit server
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
1. `python epistasis_submit.py <prefix> [options]`
	Use `ls results/<prefix> | wc -l` to check the current returned file numbers
1. run the following command on local machine: `scp -r <CONDOR_ADDRESS>:results <destination_directory_at_Parks_Lab>`

## Requirements
* CHTC account
* SQUID directory (CHTC uses uses this to transfer large files)

## Re-running jobs
If memory requirements are set too low, it can happen that some jobs will fail while most jobs will complete successfully. In these cases, it is ideal to request more memory and re-run only the failed jobs rather than re-running all the jobs.

### Procedure
1. Get the job numbers for the jobs that need to be re-run.
1. Using either [*writeFailedJobs.py*](https://github.com/Parks-Laboratory/condor_tools) or manually, create a file in the **data/** directory with one job number per line:

	```
	8605
	8689
	8707
	8805
	8810
	8836
	8840
	```

1. Call **epistasis_submit.py** with **--rerun _file_with_jobs_to_rerun_**, and whatever other flags are necessary to ensure these jobs succeed this time.

### Regarding naming of files from rerun jobs
The .out and .err files stored under condor_out will have a different cluster number than the original run, and their job/process numbers will start again at 0, but the .gwas files will be labeled with the appropriate job numbers from the file specified by the --rerun flag. Therefore, if all re-run jobs complete successfully, these .gwas files can be put in the same directory as the .gwas files from the original run, and the directory will contain a complete set of results with consistent file names.
