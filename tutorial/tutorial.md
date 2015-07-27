
# Pipeline Tutorial
For this example we will run a small data set through using a command line approach.  We are working with one MFE sample and the a wsn33 plasmid control. For computational speed these fastq files have been reduced by 75%.  We will later show how this approach can be implemented using pbs on a computing cluster.

For consistency all these commands can be run from the variant_pipeline/tutorial directory.

### 0) Dependencies getting started

	Begin by copying the variant_pipeline directory either to you scratch folder on flux or a directory on your machine.  There are up to data copies in the Box Sync folder or in /scratch/alauring_fluxm/common. This will give you your own working copy. (We can talk about git and tracking changes later)
	*Note: the pipeline has problems using paths that contain spaces. So please don't use spaces in your directory names.*

	The variant pipeline is essentially R, bowtie2, samtools, python, and picard wrapped in a bpipe pipeline. Picard and bpipe come prepacked in the libs directory so there is no need to install those.
	Bowtie2, samtools, and the R and python libraries can be installed using the instructions in the README doc.  

	If running on flux you only need to install the R packages.  The other dependencies are already available and are simply loaded as modules. (see bin/variant_pipeline.pbs for an example.) To load R on flux type

	```bash
	module load med
	module load R/3.1.1
	R
	```
	To see if you already have the needed packages you can run

	```R
	library("library_of_interest")
	```
	in R

	Ok, now that everything is set up let's get down to business.
### 1) fastq setup


	Often the fastq files we are working with with be gzipped.  If this is the case they will end in .gz and we can unzip them using the command


 ```bash
 gunzip *.gz
 ```
 This will unzip all the gzipped files in the current directory. It might take awhile if there are a lot (minutes).

 The next step will be to name the fastq file properly. Bpipe requires the fastq file to be named in the following format *sample_name.#.read_direction.fastq*, where # is the number of fastq file for the given sample and read direction (usually 1 for miseq) and read_direction is a 1 or 2 and indicates forward or reverse reads.

  Don't fret, you don't have to rename your samples by hand. To do this we'll use the change_names_miseq.py script (when working with hiseq runs naming is slightly different so we'll use the change_names_hiseq.py script). Before running the script we will test it to make sure we are naming things as we expect.  Note the default is to copy the original files to a new file. This leaves the original unchanged. Additionally, the script will not copy or move anything unless you run it with the -run flag.  Omitting this flag runs the program in test mode.  It will print what it proposes to do and make a mock log.  This ensures you don't do anything hastily.  For more information about the script simply type

 ```bash
python ../scripts/change_names_miseq.py -h
 ```
 Let's run the test
 ```bash
python ../scripts/change_names_miseq.py -s data/fastq_original/ -f data/fastq/
 ```

 Everything looks good so lets do it for real by adding the -run option

 Let's run the test
 ```bash
python ../scripts/change_names_miseq.py -s data/fastq_original/ -f data/fastq/ -run
 ```
*Note: a log of the name changes was made in fastq/renaming_log.txt for posterity*

Now that we have our samples prepped we can run the pipeline.

### 2) Running the pipeline

The pipeline is run by in [bpipe](https://code.google.com/p/bpipe/wiki/Overview) using a python wrapper.  To see our options type.

```bash
python ../bin/variantPipeline.py -h
```
We need to provide the directory containing the fastq files, the directory were we want the output  to go ( it will be made if it doesn't exist), our reference for bowtie (made above), and the name of our control sample.  To make sure everything is in working order we can run the pipeline in test mode by activating the -t option.

```bash
python ../bin/variantPipeline.py -i data/fastq/ -o worked_data -r data/reference/wsn33_wt_plasmid -p Plasmid_control -t
```

It seems like everything is in order so we'll let it rip.  This took about 5 min on my old macbook pro.

```bash
python ../bin/variantPipeline.py -i data/fastq/ -o worked_data -r data/reference/wsn33_wt_plasmid -p Plasmid_control
```


__If you are running on the flux__ Then instead of running the above command from the commad line we will run this command in a pbs script.  An example of this can be found in bin/variant_pipeline.pbs.  For larger sample sets you'll need to adjusted the memory and walltime limits.  We can use a total of 2 processors and 48 gb of mem.  A detailed description of pbs scripts can be found [here](http://arc-ts.umich.edu/software/torque/).

submit the job using

```bash
qsub ../bin/variant_pipeline.pbs
```

and check the status of it using

```bash
qstat -u yourusername
```
*You may find yourself in the unfortunate position of seeing a failed pipeline once in a while. If a certain stage fails, I have found it helpful to delete the your_output_dir/.bpipe directory, which is a hidden directory, and resubmitting the pipeline. Bpipe will rerun the pipeline from the begining, but it will skip steps whose output already exists.*


The output from a successfull run can be found in your_output_dir/mapq.all.sum.csv.
This contains the called variants and the data needed to filter them to your hearts desire.


### 3) Analyzing data

Bpipe keeps a log of all the commands it runs in 'commadlog.txt'. This can be useful for debugging.  

The pipeline does not carry out any secondary analysis. It only provides putative variants and information regarding how trustworthy those calls are.  It is up to you to sort through the putative variants (found in mapq/all.sum.csv) using Excel (booo!!) or R (Hooray!!).  Currently we are setting a p.value cut off of 0.01 (p.val < 0.01) and for most applications a frequency cut off of 0.5% (freq.var>0.5%).  Additionally we require variants to be uniformly distributed across the reads on which they are found.  We achieve this by requiring the average read position to be in the middle 50% of the read length. (Read_pos>50 & Read_pos<200, for a Miseq run with 2X250 reads).  

An example of how to begin analysis on these samples can be found in the results directory.
