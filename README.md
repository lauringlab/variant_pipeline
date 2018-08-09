# Lauring Variant Pipeline

Nominates candidate variants by comparing the sequences in a test sample to those found in a plasmid control.
The Pipeline runs as one phase which takes in fastq files and outputs putative variants as well as all base call above a set frequency. It is then up to the user to filter the putative variants based on the characteristics provided.

## Directory list

-   bin
    _ variantPipeline.py : a python wrapper that runs a provided bpipe pipeline
    _ variant_pipeline.pbs : an example pbs script used to implement the Pipeline
-   doc
    -   workflow diagram, examples
-   lib
    -   supporting libraries - picard and bpipe live here
-   packrat
    -   The R dependencies are listed in the lock file. Also the packages will be downloaded here on set up
-   scripts
    _ supporting scripts (bpipe, python, R)
    _ some of these are old and not used others are helpful for downstream analysis but are not used by the current stages
    -   most have been written to provide a useage message when called with the -h flag
-   test

    -   automated tests (mostly python testing the python pipeline)

-   tutorial
    -   The directories and instructions needed to run the tutorial. Instructions can be found in the tutorial readme

## bin/variantPipeline.py

This script is a thin python wrapper that takes in a bpipe pipeline, input files, output directory and an options yaml. Whenever this is launched, the bpipe scripts are copied from the scripts directory and stored in the output directory as a log of what was run. the output directory will be made if it doesn't exist.

Usage: python variantPipeline.py -h

    See the tutorial for more information.

    *NOTE: Your fasta is used in the variant calling step and needs to end in .fa*

## Outputs

There are 3 main pipelines that can be run. All of the stages for the pipelines are held in ./scripts/variantPipeline.bpipe.stages.groovy

### Basic alinging scripts/aligning_pipeline.groovy

-   cutadapt
    -   the trimmed fastq files - these are trimmed based on NEBnext primers which is hard coded in the stage
-   fastqc
    -   fastqc data on samples
-   align
    -   The aligned bam and sorted sam files
-   removed_duplicated
    -   bam files with duplicate reads removed

### DeepSNV pipeline scripts/deepsnv_pipeline.groovy

Runing this pipeline after the one above is the same as the old single pipeline.

-   deepSNV
    -   csv summary files, coverage files and fasta files from deepSNV
-   parsed_fa
    -   deepSNV outputs a concatenated fasta file. The parsed ones are here.
-   Variants
    -   csv files containing all variants and additional qualty data about each one. (Mapq, phred, read position ect.)
-   Filter Variants
    -   csv files containing variants that meet quality thresholds
-   Final Variants
    -   csv files containing variants that meet quality thresholds including amino acid information

### python pipeline to call all variants and sequencing errors scripts/python_pipeline.groovy

-   consensus
    -   The consesus seqeunce of each sample
-   position-stats
    -   JSON files with all bases called at every position including amino acid designation

## Dependencies

Note : _Flux is the name of the computing core used by our lab at the Univeristy of Michigan. Some of the directions may be specific to those working on this platform_

The pipeline comes with many of the required programs (bpipe and pycard); however, bowtie2, samtools and certain R and python libraries are required by the variant calling.

To run these all pipelines you must have the java developer kit installed. It can be installed from [here](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html). If bpipe doesn't run this is the first place to start.

All the other depedencies, except R and the R packages, are handled by conda. Install conda by following the tutorial [here](https://conda.io/docs/user-guide/overview.html).

We can install the conda environment with the following command (run from the variant_pipeline/ directory)

```
conda env create -f scripts/environment.yml
```

We have to activate the environment before running the commands below.

```
conda activate variant-pipeline
```

On flux we can achieve an equivalent environment by loading the following modules

```
module load muscle
module load bowtie2
module load python-anaconda2/201704
module load fastqc
module load R
```

The R modules are managed by packrat. I am using R 3.5.0. From the main directory run

```
R
packrat::restore()
```

to download the needed dependencies. They should be placed the packrat/lib directory. This is important since the R script will look for them there. You may need to install packrat first if you don't have it.

Adapted and developed by JT McCrone based on work done by
Chris Gates/Peter Ulintz
UM BCRCF Bioinformatics Core
