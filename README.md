# Lauring Variant Pipeline

Nominates candidate variants by comparing the sequences in a test sample to those found in a plasmid control.
The Pipeline runs as one phase which takes in fastq files and outputs putative variants.  It is then up to the user to filter the putative variants based on the characteristics provided (p value, frequency, and read distribution are the most useful at this point)

## Directory list
* bin
	* variantPipeline.sh : align reads, sort, call variants using deepSNV, and characterize putative variants.
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

## Workflow summary
![variantPipeline workflow image](https://github.com/umich-brcf-bioinf/lauring-variant-pipeline/blob/develop/doc/variantPipeline.png)


## bin/variantPipeline.sh (phase 1)
 This script is a thin bash wrapper around a bpipe pipeline which in turn calls fastqc, pydmx-al, bowtie, picard. Whenever this is launched, the bpipe scripts are overwrittem from the scripts directory


Usage: variantPipeline.sh {input_dir} {output_dir} {reference} {plasmid control name}

* Inputs:  
	* dir containing left fastq, right fastq named as sample.#.read_direction.fastq
		* the python scripts change_names* can be used as a fast means of renaming fastq to this format.
	* path to the reference genome for alignment
			made using
	```bash
	bowtie2-build WSN33.fa WSN33
	```
	Where WSN33.fa is your fasta file
	* plasmid control : the sample name of the plasmid control fastq.


* Outputs:
	* 01_fastqc : zips containing a brief summary report on quality of input fastqs (in development)
	* 03_align : bam files from bowtie2 alignment
	* 04_remove_duplicates : bam files with duplicates marked by picard
	* 05_Coverage : pileup and .cov (text) files containing the coverage at each basepair
	* variants : csv files of the deepSNV variant calls (bonferroni p<0.1) and .fa consensus sequences
		* an empty csv is made for the plasmid control in order to appease bpipe.
	* mapq : sum.csv files made by adding mapping quality, read position, and phred information to the putative variants found in the variant output, and reads.csv containing similar information for each read responsible for a variant
	* mapq and 05_Coverage contain files begining in "all" these files contain the output for each sample concatenated into one file for your viewing and analysing pleasure.


## Dependencies

The pipeline comes with many of the required programs (bpipe and pycard); however, bowtie2, samtools and certain R  and python libraries are required by the variant calling.

Note: *The R packages must be installed under your username on Flux.  The other dependencies are simply added as modules.*

* R (installation instructions to be run in R)
 	* deepSNV

		```R
		source("http://bioconductor.org/biocLite.R")
		biocLite("deepSNV")
		```
	* plyr

		```R
		install.packages("plyr")
		```
	* reshape2

	```R
	install.packages("reshape2")
	```

* python	(installed through pip in terminal)
	* pysam 0.8.2.1

	```bash
	pip install pysam
	```
* bowtie2
	[Getting started with bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)


## Pipeline Example
For this example we will run a small data set through using a command line approach.  We will later show how this approach can be implemented using pbs on a computing cluster.





Adapted and developed by JT McCrone based on work done by
Chris Gates/Peter Ulintz
UM BCRCF Bioinformatics Core
