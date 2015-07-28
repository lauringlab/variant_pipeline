
# Pipeline Tutorial
For this example we will run a small data set through using a command line approach.  We are working with one MFE sample and the a wsn33 plasmid control. For computational speed these fastq files have been reduced by 75%.  This approach can be run on flux as is, or run through a pbs for larger data sets.

For consistency all these commands can be run from the variant_pipeline/tutorial directory.

### 0) Getting started

Log onto the flux platform  by typing the following in terminal.

```
ssh your_username@flux-login.engin.umich.edu
```

 You will be prompted for your MToken.  No characters will appear as you type.  You make a mistake you will be locked out until you bring cookies to lab.  Just kidding, but you do have to wait until your MToken refreshes to a new number.

 You will then be asked for your level one.  Again no characters appear as you type.

 Once on flux you will automatic all begin in your home directory (~/ which is a shortcut for /home2/username/).  We have limmited space in these directories so we will typically work from scratch directories which provide more memory for active work but should not be used for longterm storage.  We will add the variant pipeline to our home directory so it is easily accessible from anywhere.

 You are reading this tutorial so you must be on github.  Our first task will be to clone the github repository (the variant caller in all its glory) to you home directory on flux.
 
 On flux type
 
 ```
 git clone https://username:password@github.com/jtmccr1/variant_pipeline.git
 ```
where username:password are your github username and password.
 
 Now we will make sure all scripts are exicutable.  

 ```
 ls -l
 ```

The output should look somthing like this

```
drwx------ 9 mccrone microbio     4096 Jul 27 13:44 variant_pipeline
```

The first characters mean this is directory (d) and the owner (mccrone) has read (r), write (w), and execute (x) permissions. If you don't have an x in these first 4 letters you can add it with

 ```
 chmod -R +x variant_pipeline
 ```

 (Again the -R is a recursive option) You need to do this to be able to run the scripts that come with the pipeline.

Now we are ready to begin the tutorial. Let's go there now.

```
cd variant_pipeline/tutorial
```

The variant pipeline is essentially R, bowtie2, samtools, python, and picard wrapped in a bpipe pipeline. Picard and bpipe come prepacked in the libs directory so there is no need to install those. Samtools and bowtie2 are already on flux.  So all we need to do is install the deepSNV package per the instructions in the README doc.

To see if you already have deepSNV you can run
```
module load med
module load R/3.1.1
R
library("deepSNV")
```
If you don't have it you can install it using

```
source("http://bioconductor.org/biocLite.R")
biocLite("deepSNV")
```
while still in R.  It should take a while to install and you will be prompted to install in your private packages.  Please respond y when prompted.

use 

```
q()
```
to exit R when you are done.

Next, we will need to load the needed programs using


```
module add med
module add fastqc
module add java
module add bowtie2
module add samtools
module add python
module add R/3.1.1
module add pysam/0.8.2.1
```

Ok, now that everything is set up let's get down to business.
### 1) fastq setup



Ok, now that everything is set up let's get down to business.
### 1) fastq setup



```bash
 gunzip *.gz
 ```
This will unzip all the gzipped files in the current directory. It might take awhile if there are a lot (minutes).

The next step will be to name the fastq file properly. Bpipe requires the fastq file to be named in the following format *sample_name.#.read_direction.fastq*, where # is the number of fastq file for the given sample and read direction (usually 1 for miseq) and read_direction is a 1 or 2 and indicates forward or reverse reads.

Don't fret, you don't have to rename your samples by hand. To do this we'll use the change_names_miseq.py script (when working with hiseq runs naming is slightly different so we'll use the change_names_hiseq.py script). Before running the script we will test it to make sure we are naming things as we expect.  Note the default is to copy the original files to a new file. This leaves the original unchanged. Additionally, the script will not copy or move anything unless you run it with the -run flag.  Omitting this flag runs the program in test mode.  It will print what it proposes to do and make a mock log.  This ensures you don't do anything hastily.  For more information about the script simply type

```
python ../scripts/change_names_miseq.py -h
```
Let's run the test

*Note: if the final directory ("data/fastq" in this case) doesn't exist it will be made*
```
python ../scripts/change_names_miseq.py -s data/fastq_original/ -f data/fastq/
```
Everything looks good so lets do it for real by adding the -run option

```
python ../scripts/change_names_miseq.py -s data/fastq_original/ -f data/fastq/ -run
```
*Note: a log of the name changes was made in fastq/renaming_log.txt for posterity*

Now that we have our samples prepped we can run the pipeline.

### 2) Running the pipeline

The pipeline is run by in [bpipe](https://code.google.com/p/bpipe/wiki/Overview) using a python wrapper.  To see our options type.

```
python ../bin/variantPipeline.py -h
```

We need to provide the directory containing the fastq files, the directory were we want the output  to go ( it will be made if it doesn't exist), our reference for bowtie (provided by the tutorial, see README on how to make one for your sample set), and the name of our control sample.  To make sure everything is in working order we can run the pipeline in test mode by activating the -t option.

```
python ../bin/variantPipeline.py -i data/fastq/ -o worked_data -r data/reference/wsn33_wt_plasmid -p Plasmid_control -t
```

It seems like everything is in order so we'll let it rip.  This took about 5 min on my old macbook pro. On the flux it won't take long at all.

```
python ../bin/variantPipeline.py -i data/fastq/ -o worked_data/ -r data/reference/wsn33_wt_plasmid -p Plasmid_control
```


__If you are running on the flux__ Then instead of running the above command from the command line we would usually run this command in a pbs script. This is because running memory intensive commands on the login node slows everyone down.  A pbs script tells flux to set aside a separate node just for our work. An example of this can be found in bin/variant_pipeline.pbs.  For larger sample sets you'll need to adjusted the memory and walltime limits.  We can use a total of 2 processors and 48 gb of mem.  A detailed description of pbs scripts can be found [here](http://arc-ts.umich.edu/software/torque/).
The script can be edited using nano.

Let's run the same pipeline using a pbs script.

The first lines of a pbs script begin with "#PBS" and give the scheduler information regarding our job.
Let's make some modifications using nano.

```
nano ../bin/variant_pipeline.pbs
```

*Note you can't use the mouse to navigate and you can save by pressing cntrl+x*

The line

> #PBS -N test_PR8

signifies the name of the job.  Let's change "test_PR8" to "tutorial"

> #PBS -M mccrone@umich.edu

Tells the computer to send you an email at the start and end of the job let's change it so I don't get all your emails.

> #PBS -l nodes=1:ppn=2,mem=24gb,walltime=2:00:00

signifies how many nodes, processors per node, memory, and max time we want the job to run.  For this small job lets use 10gb of memory and limit the wall time to 10 min 00:10:00.


We then can submit the job using

```bash
qsub ../bin/variant_pipeline.pbs
```


and check the status of it using

```bash
qstat -u yourusername
```

*You may find yourself in the unfortunate position of seeing a failed pipeline once in a while. If a certain stage fails, I have found it helpful to delete the your_output_dir/.bpipe directory, which is a hidden directory, and resubmitting the pipeline. Bpipe will rerun the pipeline from the begining, but it will skip steps whose output already exists.*


A log of the job output can be found at tutorial.o########## where tutorial is the name of the job and ####### is the job Id.  We can page through the output using

```
more tutorial.o*
```
and spacebar to page down.

At the bottom we should find

>---------------------  Pipeline Success  -------------------------

The output data from a successful run can be found in your_output_dir/mapq.all.sum.csv.
This contains the called variants and the data needed to filter them to your hearts desire.

Additionally Bpipe keeps a log of all the commands it runs in 'commandlog.txt'. This can be useful for debugging.  

### 3) Analyzing data



The pipeline does not carry out any secondary analysis. It only provides putative variants and information regarding how trustworthy those calls are.  It is up to you to sort through the putative variants (found in mapq/all.sum.csv) using Excel (booo!! :-1:) or R (Hooray!! :bowtie:).  Currently we are setting a p.value cut off of 0.01 (p.val < 0.01) and for most applications a frequency cut off of 0.5% (freq.var>0.5%).  Additionally we require variants to be uniformly distributed across the reads on which they are found.  We achieve this by requiring the average read position to be in the middle 50% of the read length. (Read\_pos>50 & Read\_pos<200, for a Miseq run with 2X250 reads).  

To analyze the data transfer the mapq/all.mapq.* and 05_Coverage/all.cov.csv to your own machine using your favorite file transfer software.  Cyberduck and Filezilla are good places to start.

An example of how to subset the data and plot coverage can be found in the results directory of the tutorial.  You can look at this my opening it from the box sync folder.


### 4) Modifying for your own data

So you just got some illumina data back! Bully for you! How to analyze it.  Using either cyberduck or better still globus (see command line tools in box synce) transfer your run data to the appropriate directory on NAS ( which is backedup). The path to our NAS is "/nfs/med-alauring" and there are directories for raw data that are organized by year. Put your raw data in the appropriate directory.

Next cd the scratch directory where we have more memory to run our large jobs

```
cd /scratch/alauring_fluxm/
ls
```
Look there is folder just for you! cd into it and we can begin.

## setup 
Let's setup an experimental directory
```
mkdir exp_label
cd exp_label
mkdir data
mkdir scripts
```

*Note you may have to make a reference file for bowtie to align to.  I like to keep mine in data/reference.  You can use the command in the readme file to make your reference so long as you already have a fasta file.*


Now we'll rename and move our fastq's from the NAS to our data directory.  We just have to tell the computer where to find our scripts.  These commands should look familar.

```
python ~/variant_pipeline/scripts/change_names_miseq.py -s path/to/data/on/NAS -f data/fastq/ 
```
If this looks good let's run it

```
python ~/variant_pipeline/scripts/change_names_miseq.py -s path/to/data/on/NAS -f data/fastq/ -run
```

Now we can copy the pbs script from the tutorial and modify it to suit our purposes.  

```
cp ~/variant_pipeline/bin/variant_pipeline.pbs
```

The last line should read

```
python ~/variant_pipeline/bin/variantPipeline.py -i data/fastq/ -o worked_data/ -r path/to/reference/name -p your_plasmid_control
```

Then just submit using qsub as before and sit back while the computer does the rest. :smiley:

