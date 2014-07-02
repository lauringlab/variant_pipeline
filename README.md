# Lauring Variant Pipeline

Nominates candidate variants by comparing pairs of flu strain sequences.
The pipeline runs in two phases:

* a batch phase accepts paired fastq sequences specially prepared for pydmx deduplication (see details on pydmx below). This phase deduplicates/demultiplexes individual barcoded samples emitted a set of alignmed bams. 
* a "manual" calling phase accepts a test and control bams and calls the variants sing deepSNV

## Directory list
* bin
	* variantPipeline.sh : reduce/demultiplex all barcoded samples from a fastq pair
	* deepSNV.sh : call variants from individual test and control pair
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
Accepts a directory containing a barcode file (e.g. bars.csv) and pair of barcoded fastqs and emits a demultiplexed set of fastq pairs representing the consensus sequence. This script is a thin bash wrapper around a bpipe pipeline which in turn calls fastqc, pydmx-al, bowtie, picard. Whenever this is launched, the bpipe scripts are overwrittem from the scripts directory


Usage: variantPipeline.sh {input_dir} {output_dir}

* Inputs:  
	* dir containing left fastq, right fastq and bars.csv
	* reference genome [hardcoded in scripts/variantPipeline.bpipe.config]

* Outputs:
	* 01-fastqc : zips containing a brief summary report on quality of input fastqs
	* 02-pydmx : several directories that hold the intermediate and final results of pydmx processing; see details on pydmx-al below
	* 03-align : bam from bowtie2 alignment
	* 04-mark_duplicates : bam files with duplicates marked by picard


## bin/deepSNV.sh (phase 2)
Accepts a reference genome, test bam, and a control bam and emits variants called by deepSNV in CSV and VCF formats. This is a thin bash wrapper around the deepSNV.r script.

Usage: deepSNV.sh {reference_fasta} {test.bam} {control.bam} {output_dir}

## pydmx-al details
Deduplicates PCR replicates. This extends the original pydmx implementation by adding steps to trim reads of unequal length and demultiplex the reads based on barcode sequence.

### details on deduplication 
First filters to the set of fastq sequences where specific positions matched known barcode and then (based on PCR-counter tag) collapses duplicates into a single sequence.
* Inputs
 * forward and reverse paired-end fastq files
 * csv of expected barcodes
* Outputs 
 * the subset of deduplicated sequences, trimming off pcr-counter tag and sample barcode, and replacing the sequence header with the barcode and duplication count.
 * a text summary file showing barcodes, duplication, and sequences

Looking for barcodes that include **CACTGT**, given *1.fastq*:
> @M02127:15:000000000-A6WCU:1:1101:16578:1599 1:N:0:1
> TGGAA<b>CACTGT</b>TGGCCCTGTCCATTTTAGAAACCAAGTCAAAATACGTCGG...<br/>
> \+ <br/>
> D1111B1>D1FD111AAF1FBF3DAFF30133BB000BA2211B1B0BE//...


pydmx would emit *deduped.1.fastq*:
> @<b>CACTGT</b>:DUPE_3:ID_1047:FLAG_1 1
> TGGCCCTGTCCATTTTAGAAACCAAGTCAAAATACGTCGGAGAGTTGACA...<br/>
> \+ <br/>
> D111AAF1FBF3DAFF30133BB000BA2211B1B0BE//////BD12B1...

As well as the text equivalent *deduped.summary.txt* where each row represents a deduplicated, matched pair: 

pydmx_id | barcode | left_key | right_key | num_dupes | barcode_flag | left_seq | left_qual | right_seq | right_qual | duplicate_ids
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
1047 | CACTGT | <b>TGGAACACTGT</b>TGGCCCT | <b>CAACACTTCAA</b>TGTACAC | 3 | 1 | TGGCCCTGTC... | ... | TGTACACACTGCT... | ... | [M02127:15:000000000-A6WCU:1:1101:16578:1599]...


Note:
* The original PCR-counter (TGGAA) and sample barcode (CACTGT) have been removed in the output sequence.
* The output sequence identifier "@CACTGT:DUPE_3:ID_1047:FLAG_1 1" represents:
 * CACTGT : The matching barcode
 * DUPE_3 : Indicates this sequence appeared three times in the original file
 * ID_1047 : An arbitrary id which matches between left, right, and summary files.
 * FLAG_1 : Indicates a barcode was found in the forward (left) sequence. FLAG_2 indicates barcode in right strand and FLAG_3 in both strands.
 * 1 : This is the forward (left) fastq (2 = reverse strand, 3 = both)
* If a barcode is found on one side but not the other, pydmx will look for and remove the sequencing primer on the side missing the barcode. (This affects a very small percentage of reads which the upstream trimmer missed.)
* Pydmx cuts the paired input files into chunks to parallelize across multiple processors. Because it chunks the files by bytes and not lines, the input files *must* be exactly the same size. Mismatched files are rejected immediately.


See pydmx source for more information on the expected structure of the PCR-counter and sample barcode tags.

---

Chris Gates/Peter Ulintz
UM BCRCF Bioinformatics Core
