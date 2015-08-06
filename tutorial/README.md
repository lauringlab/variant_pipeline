
# Pipeline Tutorial
In this example we will run a small data set through using a command line approach. We are working with one MFE sample and a wsn33 plasmid control. For computational speed, these fastq files have been reduced by 75%. This approach can be run on flux as is, or run through a pbs for larger data sets. The last section provides instructions for working with real data.

For consistency all of the commands can be run from the variant_pipeline/tutorial directory. If you are running this from a different directory (for example, when you are working with your real data on scratch in section 4), make sure the paths to your files are correct. You may need to edit them in the commands provided.

## Table of Contents

* [Before you begin](#before-you-begin)
* [0) Getting started](#getting-started)
* [1) Fastq setup](#fastq-setup)
* [2) Running the pipeline](#running-the-pipeline)
* [3) Analysis](#analysis)
* [4) Working with real data](#working-with-real-data)

<a name="before-you-begin"/>
## Before you begin
You will need:
- Access to the lab's allocation on flux. Adam will need to email the administrator at hpc.
- An MToken for two factor authentication. See [MSIS](http://www.mais.umich.edu/mtoken/mtoken_distribution.html).
- Know basic unix commands. See tutorials and lists on Mbox/Lauring Lab/Command Line Tools.
- Know the basics of flux organization and access. See Ten easy steps in MBox/Lauring Lab/Command Line Tools.



<a name="getting-started"/>
## 0) Getting started

Log onto the flux platform  by typing the following in a terminal window.

```
ssh your_username@flux-login.engin.umich.edu
```

 You will be prompted for your MToken. No characters will appear as you type. If you make a mistake, you will be locked out until you bring cookies to lab. Just kidding, but you do have to wait until your MToken refreshes to a new number.

 You will then be asked for your level one. Again no characters appear as you type.

 Once on flux you will automatically begin in your home directory (~/ which is a shortcut for /home/username/ ). If you want to check your location, simply type the unix command "pwd" and you should see /home/your_username. We have limited space in these directories (80gb) so we will typically work from scratch directories which provide more memory for active work. However, scratch should not be used for longterm storage. Therefore, we will add the variant pipeline to our home directory so it is easily accessible from anywhere.

 You are reading this tutorial so you must be on github. Our first task will be to clone the github repository (the latest version of the variant caller in all its glory) to you home directory on flux.
 
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

(The -R is a recursive option) You need to do this to be able to run the scripts that come with the pipeline.

Although the pipeline is fully functional, we may make slight modifications in the future. To ensure you are working with the most up to date version you can use 

```
git pull
```

This should be executed from somewhere in  the variant_pipeline directory and  will "pull" the most update version from github.  If you have modified any of the files (and this is unlikely) you will be asked to commit those changes before you can execute "git pull". Without getting too much into git.

```
git commit -am "committing before a pull"
```

will commit your changes and label them "committing before a pull". You can now pull.

Now the code for the variant pipeline is installed in your home directory in a sub-directory called variant_pipeline. We are ready to begin the tutorial. Let's go there now, by moving to a subdirectoy within the variant_pipeline directory.

```
cd variant_pipeline/tutorial
```

The variant pipeline is essentially R, bowtie2, samtools, python, and picard wrapped in a bpipe pipeline. Picard and bpipe come prepacked in the libs directory so there is no need to install those. Samtools and bowtie2 are already on flux. So all we need to do is install the deepSNV package per the instructions in the README doc.

To see if you already have deepSNV you can run
```
module load med
module load R/3.1.1
R
library("deepSNV")
```
If you are in R, you should see the cursor change to ">". If you get a message saying you don't have deepSNV  you can install it using

```
source("http://bioconductor.org/biocLite.R")
biocLite("deepSNV")
```
while still in R.  It should take a while to install and you will be prompted to install in your private packages. Please respond y when prompted. You may also be prompted to update all/some/none of the packages. Respond a when prompted.

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

Ok, now that everything is set up, let's get down to business.

<a name="fastq-setup"/>
## 1) Fastq setup

The first step is to name the fastq file properly. They arrive from the sequencing core with names that Bpipe can't make sense of. Bpipe requires the fastq file to be named in the following format *sample_name.read_direction.#.fastq*, where # is the number of fastq file for the given sample and read direction (usually 1 for miseq) and read_direction is a 1 or 2, indicating forward or reverse reads.

Don't fret, you don't have to rename your samples by hand. To do this we'll use the change_names_miseq.py script (when working with hiseq runs naming is slightly different so we'll use the change_names_hiseq.py script). Before running the script, we will test it to make sure we are naming things as we expect. Note that the script will copy the original files to a new file. This leaves the original unchanged. Additionally, the script will not copy anything unless you run it with the -run flag.  Omitting this flag runs the program in test mode. It will print what it proposes to do and make a mock log. This ensures you don't do anything hastily. For more information about the script simply type
 		 
 ```		
 python ../scripts/change_names_miseq.py -h		 
 ```

 		 
*Note: if the final directory ("data/fastq" in this case) doesn't exist it will be made. __Also, fastq files are gzipped when we get them from the sequencing core (they end in .gz). This script will detect files that end in ".fastq" and ".fastq.gz".  It will copy the unzipped and gzipped files from the -s directory (data/fastq_original/)  to  -f (data/fastq/) and then it will unzip all zipped files in -f so that we can use them in analysis.  This may take some time for large files__*

Let's run the test.

```
python ../scripts/change_names_miseq.py -s data/fastq_original/ -f data/fastq/
```
Everything looks good so lets do it for real by adding the -run option

```
python ../scripts/change_names_miseq.py -s data/fastq_original/ -f data/fastq/ -run
```
*Note: a log of the name changes was made in fastq/renaming_log.txt for posterity*



Now that we have our samples prepped we can run the pipeline.


<a name="running-the-pipeline"/>
## 2) Running the pipeline

The pipeline is run by in [bpipe](https://code.google.com/p/bpipe/wiki/Overview) using a python wrapper.  To see our options type.

```
python ../bin/variantPipeline.py -h
```

We need to provide the directory containing the fastq files, the directory where we want the output to go (it will be made if it doesn't exist), our reference for bowtie (provided by the tutorial, see README on how to make one for your sample set), and the name of our control sample. To make sure everything is in working order we can run the pipeline in test mode by activating the -t option.

```
python ../bin/variantPipeline.py -i data/fastq/ -o worked_data -r data/reference/wsn33_wt_plasmid -p Plasmid_control -t
```

It seems like everything is in order so we'll let it rip. This took about 5 min on my old macbook pro. On the flux it won't take long at all.

```
python ../bin/variantPipeline.py -i data/fastq/ -o worked_data/ -r data/reference/wsn33_wt_plasmid -p Plasmid_control
```


__If you are running on the flux__ Then instead of running the above command from the command line we would usually run this command in a pbs script. This is because running memory intensive commands on the login node slows everyone down.  A pbs script tells flux to set aside a separate node just for our work. An example of a pbs script can be found in bin/variant_pipeline.pbs.  For larger sample sets you'll need to adjust the memory and walltime limits___can you provide suggested times and memory for a hiseq run or a miseq run?___.  We can use a total of 2 processors and 48 gb of mem.  A detailed description of pbs scripts can be found [here](http://arc-ts.umich.edu/software/torque/).
The script can be edited using the easy text editor, nano.

Let's run the same pipeline using a pbs script.

The first lines of a pbs script begin with "#PBS" and give the scheduler information regarding our job.
Let's make some modifications using nano.

```
nano ../bin/variant_pipeline.pbs
```

*Note you can't use the mouse to navigate in nano, but you can use the arrow keys. You can save by pressing cntrl+x*

The line

> #PBS -N test_PR8

signifies the name of the job. Let's change "test_PR8" to "tutorial"

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

<a name="analysis"/>
## 3) Analysis



The pipeline does not carry out any secondary analysis. It only provides putative variants and information regarding how trustworthy those calls are.  It is up to you to sort through the putative variants (found in mapq/all.sum.csv) using Excel (booo!! :-1:) or R (Hooray!! :bowtie:).  Currently we are setting a p.value cut off of 0.01 (p.val < 0.01) and for most applications a frequency cut off of 0.5% (freq.var>0.5%).  Additionally we require variants to be uniformly distributed across the reads on which they are found.  We achieve this by requiring the average read position to be in the middle 50% of the read length. (Read\_pos>50 & Read\_pos<200, for a Miseq run with 2X250 reads).  

To analyze the data transfer the mapq/all.mapq.* and 05_Coverage/all.cov.csv to your own machine using your favorite file transfer software. Cyberduck and Filezilla are good places to start.

An example of how to subset the data and plot coverage can be found in the results directory of the tutorial. You can look at this by opening it from the box sync folder.

<a name="working-with-real-data"/>
## 4) Working with real data

So you just got some Illumina data back! Bully for you! Now to analyze it. Using either cyberduck, or better still, [globus](http://arc.umich.edu/flux-and-other-hpc-resources/flux/using-flux/transferring-files-with-globus-gridftp/) transfer your run data to the appropriate directory on NAS (which is backed up regularly by the University and also backed up to a lab external hard drive). DO NOT delete your data from the portable hard drive unless you check with Adam first. NEVER edit these primary sequence files. The path to our NAS is "/nfs/med-alauring" and there are directories for raw data that are organized by year. (See the flux directory sctructure below). It is a good idea to rename your directory prior to transfer so that it does not have spaces. Stay tuned for a uniform nomenclature for our lab.

### Flux directory structure

![Dir structure](https://github.com/jtmccr1/variant_pipeline/blob/master/doc/Flux.jpeg)

Once your data is in NAS, be a good neighbor and make the data accessible to everyone in the lab. The first command makes alauring-lab the group for the files.  We have to do this because the default for some people in the lab is internalmed and for others is micro.  This makes the group something everyone belongs to.  The next command gives those in the group read, write,  and execute permission.

```
chgrp -R alauring-lab path/to/your_data
chmod -R g=rwx
```

Now that the data in up on NAS. Let's get a directory set up on /scratch/ where we will do our actual work.

```
cd /scratch/alauring_fluxm/
ls
```
Look there is folder just for you! cd into it and we can begin.

### 4.1 Setup 

Let's setup an experimental directory. This will hold all of the files and data that you use for a given experiment or flux run. From your scratch directory run, where "exp_label" is the name you choose for this experiment. Make it something that provides information about the experiment and/or a date.
```
mkdir exp_label
cd exp_label
mkdir data
mkdir data/fastq
mkdir data/reference
mkdir scripts
```
You can now navigate through the exp_label directory and sub-directories to see that there is a directory called "data" that contains sub-directories for fastq files and reference files. There is also a sub-directory called "scripts" that you will rarely access.

*Note you may have to make a reference file for bowtie to align to.  I like to keep mine in data/reference.  You can use the command in the readme file to make your reference so long as you already have a fasta file. __IT MUST END IN .fa FOR THE VARIANT CALLER TO RECOGNIZE IT__*

### 4.2 Running
Now we'll rename and copy our fastq's from the NAS to our data directory.  We just have to tell the computer where to find our scripts.  These commands should look familar, and can be run from your experiment folder on scratch.

```
python ~/variant_pipeline/scripts/change_names_miseq.py -s path/to/data/on/NAS -f data/fastq/ 
```
The only differences between this and what we did above are the paths to the scripts and the data. The path to your data on NAS should be something like "/nfs/med-alauring/raw_data/2015/filename". Because the pipeline is in your home directory you can easily access it with the shortcut "~/"

If this looks good let's run

```
python ~/variant_pipeline/scripts/change_names_miseq.py -s path/to/data/on/NAS -f data/fastq/ -run
```

*Note if you data on NAS ends in .gz then it is gzipped.  This script is able to copy and unzip zipped files automatically.*


Now we can copy the pbs script from the tutorial and modify it to suit our purposes.  

```
cp ~/variant_pipeline/bin/variant_pipeline.pbs ./experiment_name.pbs
```

You can use the command "nano experiment_name.pbs" to open the text editor and edit the pbs script. Make sure you edit the path to the variantPipeline.py script by directing it to your home directory (~/). The last line should read

```
 python ~/variant_pipeline/bin/variantPipeline.py -i data/fastq/ -o worked_data/ -r path/to/reference/name -p your_plasmid_control
```

Then just submit using qsub as before and sit back while the computer does the rest. :smiley:

