# Epistasis Pipeline

## Running the Epistasis Pipeline on UW-Madison Clusters

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


NOTE:

* Sometimes the "condor_q" command fails and quits the program unexpectedly. In this case, rerun the pipeline  
* Ensure that files scripts/fastlmmc, scripts/plink, fixpheno.sh are executable (green). If not run the command `chmod +x FILE_NAME`

REQUIREMENTS:

* CHTC account
* SQUID directory

## Future work
* Modify epistasis_server.py to continuously monitor for failed jobs and to
re-submit them if appropriate
	* add periodic_release code
	* use HoldReasonCode?
* Monitor for completed jobs, delete their associated .err, .out files if appropriate to reduce space used on submit server
