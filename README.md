# Lauring Variant Pipeline

Nominates candidate variants by comparing the sequences in a test sample to those found in a plasmid control.
The Pipeline runs as one phase which takes in fastq files and outputs putative variants.  It is then up to the user to filter the putative variants based on the characteristics provided (p value, frequency, and read distribution are the most useful at this point)

## Directory list
* bin
	* variantPipeline.py : align reads, sort, call variants using deepSNV, and characterize putative variants.
	* variant_pipeline.pbs : an example pbs script used to implement the Pipeline
* doc
	* workflow diagram, examples
* lib
	* supporting libraries (mostly R libraries)
* scripts
	* supporting scripts (bpipe, python, R)
	* note that variantPipeline.bpipe.config holds several default config variables
* test
	* automated tests (mostly python)

* tutorial
		* The directories and instructions needed to run the tutorial. Instructions can be found in tutorial.html.

## Workflow summary
![variantPipeline workflow image](doc/workflow.png)


## bin/variantPipeline.py
 This script is a thin python wrapper around a bpipe pipeline which in turn calls fastqc, pydmx-al, bowtie, picard. Whenever this is launched, the bpipe scripts are overwrittem from the scripts directory and stored in the output directory as a log of what was run.


Usage: variantPipeline.py -i {input_dir} -o {output_dir} -r {reference} -p {plasmid control name}

* Inputs:  
	* -i dir containing left fastq, right fastq named as sample.read_direction[1,2].#.fastq
		* the python scripts change_names* can be used as a fast means of renaming fastq to this format.
	* -r path to the reference genome for alignment made using

	```bash
	bowtie2-build WSN33.fa WSN33
	```
	Where WSN33.fa is your fasta file
	* -p plasmid control : the sample name of the plasmid control fastq.

	See the tutorial for more information.


* Outputs:
	* __01_fastqc__ : zips containing a brief summary report on quality of input fastqs (in development)
	* __03_align__ : bam files from bowtie2 alignment
	* __04_remove_duplicates__ : bam files with duplicates marked by picard
	* __05_Coverage__ : pileup and .cov (text) files containing the coverage at each basepair
	* __variants__ : csv files of the deepSNV variant calls (bonferroni p<0.1) and .fa consensus sequences
		* an empty csv is made for the plasmid control in order to appease bpipe.
	* __mapq__ : sum.csv files made by adding mapping quality, read position, and phred information to the putative variants found in the variant output, and reads.csv containing similar information for each read responsible for a variant
	* mapq and 05_Coverage contain files begining in "all" these files contain the output for each sample concatenated into one file for your viewing and analyzing pleasure.
	* __doc__ : Bpipe makes a cool/qausi-useful html report in doc/index.html
	* __.bpipe__: a hidden directory that contains information bpipe uses to restart failed runs.  *Note: It is sometimes useful to delete this and start fresh*


## Dependencies

The pipeline comes with many of the required programs (bpipe and pycard); however, bowtie2, samtools and certain R  and python libraries are required by the variant calling.

Note: *The R package deepSNV may need to be installed under your username on Flux.  The other dependencies are simply added as modules.*

To open R on flux simply type
	```
	module load med
	module load R\3.1.1
	R
	```

	This can be done from any directory.
* deepSNV

		```
		source("http://bioconductor.org/biocLite.R")
		biocLite("deepSNV")
		```

You will be prompted to install in a local directory beginning with ~/. This means you are installing in your home directory and the library will be available just to you.  Installation should take a while.

Adapted and developed by JT McCrone based on work done by
Chris Gates/Peter Ulintz
UM BCRCF Bioinformatics Core
