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
![variantPipeline workflow image](doc/workflow.png)


## bin/variantPipeline.py
 This script is a thin bash wrapper around a bpipe pipeline which in turn calls fastqc, pydmx-al, bowtie, picard. Whenever this is launched, the bpipe scripts are overwrittem from the scripts directory


Usage: variantPipeline.sh {input_dir} {output_dir} {reference} {plasmid control name}

* Inputs:  
	* dir containing left fastq, right fastq named as sample.#.read_direction.fastq
		* the python scripts change_names* can be used as a fast means of renaming fastq to this format.
	* path to the reference genome for alignment made using

	```bash
	bowtie2-build WSN33.fa WSN33
	```
	Where WSN33.fa is your fasta file
	* plasmid control : the sample name of the plasmid control fastq.


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

Note: *The R packages may need to be installed under your username on Flux.  The other dependencies are simply added as modules.*

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
	[How to install bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)

* samtools
	[Installing samtools](http://www.htslib.org)


## Pipeline Example
For this example we will run a small data set through using a command line approach.  We are working with one MFE sample and the a wsn33 plasmid control. For computation speed these fastq files have been reduced 75%.  We will later show how this approach can be implemented using pbs on a computing cluster.

For consistency all these commands will be run from the variant_pipeline/tutorial directory.

### 1) fastq setup


	Often the fastq files we are working with with be gzipped.  If this is the case they will end in .gz and we can unzip them using the command


 ```bash
 gunzip *.gz
 ```
 This will unzip all the gzipped files in the current directory. It might take awhile if there are a lot (minutes).

 The next step will be to name the fastq file properly.  To do this we'll use the change_names_miseq.py script. Before running the script we will test to make sure we are naming things as we expect.  Note the default is to copy the original files to a new file. This leaves the original unchanged. For more information about the script simply type

 ```bash
python ../scripts/change_names_miseq.py -h
 ```
 Let's run the test
 ```bash
python ../scripts/change_names_miseq.py -s fastq_original/ -f fastq/
 ```

 Everything looks good so lets do it for real by adding the -run option

 Let's run the test
 ```bash
python ../scripts/change_names_miseq.py -s fastq_original/ -f fastq/ -run
 ```
*Note: a log of the name changes was made in fastq/renaming_log.txt for posterity*

### 2) Running the pipeline

The pipeline is run by in [bpipe](https://code.google.com/p/bpipe/wiki/Overview) using a python wrapper.  To see our options type.
```bash
python ../bin/variantPipeline.py -h
```
So we need to provide the directory containing the fastq files, the directory we were we want the output ( it will be made if it doesn't exist), our reference for bowtie (made above), and the name of our control sample.  To make sure everything is in working order we can run the pipeline in test mode by activating the -t option.

```bash
python ../bin/variantPipeline.py -i fastq/ -o worked_data -r reference/wsn33_wt_plasmid -p Plasmid_control -t
```

It seems like everything is in order so we'll let it rip.  This took about 5 min on my old macbook pro.

```bash
python ../bin/variantPipeline.py -i fastq/ -o worked_data -r reference/wsn33_wt_plasmid -p Plasmid_control
```

### 3) Analyzing data

Bpipe keeps a log of all the commands it runs in 'commadlog.txt'. This can be useful for debugging.  

The pipeline does not carry out any secondary analysis. It only provides putative variants and information regarding how trustworthy those calls are.  It is up to you to sort through the putative variants (found in mapq/all.sum.csv) using Excel (booo!!) or R (Hooray!!).  Currently we are setting a p.value cut off of 0.01 (p.val<0.01) and for most applications a frequency cut off of 0.5% (freq.var>0.5%).  Additionally we require variants to be uniformly distributed across the reads on which they are found.  We achieve this by requiring the average read position to be in the middle 50% of the read length. (Read_pos>50 & Read_pos<200, for a Miseq run with 2X250 reads).  

An example of how to begin analysis on these samples can be found in the results directory.


Adapted and developed by JT McCrone based on work done by
Chris Gates/Peter Ulintz
UM BCRCF Bioinformatics Core
