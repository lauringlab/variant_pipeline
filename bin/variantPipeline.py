import os
import os.path
import argparse
import shutil
import subprocess as s
import yaml
import sys

parser = argparse.ArgumentParser(description='This is a wrapper to set up and run the bpipe command')

parser.add_argument('options', metavar='options', nargs='+',
                    help='A YAML options file mimicing the one found in the bin directions')
parser.add_argument('-t',action='store_true',dest='test',default=False,help='Boolean switch to run program in test mode. Everything will be set up but bpipe will run in test mode')

args=parser.parse_args()

with open(args.options[0], 'r') as stream:
    try:
        options=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise "YAML error"


## Give the input arguments better names ##
input_dir=os.path.abspath(options["input_dir"])
output_dir=os.path.abspath(options["output"])
ref=os.path.abspath(options["ref"])
control=options["control"]
disp=options["disp"]
p_cut=options["p"]
method=options["method"]
open_reading=os.path.abspath(options["open_reading"])
bin_dir=os.path.dirname(os.path.realpath(__file__)) # The path to this file so we can find the scripts and lib
script_dir=os.path.abspath(bin_dir+'/..'+'/scripts/')# The path to the scripts dir relative to this location
lib_dir=os.path.abspath(bin_dir+'/..'+'/lib/') # The path to the lib dir relative to this location
bpipe_command=lib_dir+'/bpipe-0.9.8.7/bin/bpipe' # The path to the bpipe command relative the lib dir.
test=args.test
#bam=args.bam

print "Processing fastqs from " + input_dir
print "Results will be saved to " + output_dir
print "Using " + ref +" for a reference and \n" + control + " as the control sample"

git_command="git rev-list HEAD |head -n 1"
print "using version : "
s.call(git_command,shell=True)

## If the output dir does not exist make it 
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.chdir(output_dir)

# Copy the stages script and the pipeline to the output dir
shutil.copy(script_dir+'/variantPipeline.bpipe.stages.groovy',output_dir)
# if bam==True: # copy the pipeline that starts post alignment if the bam flag is set.
#     shutil.copy(script_dir+'/variantPipeline.postalign.bpipe.groovy',output_dir)
# else: # otherwise just keep keepin on
shutil.copy(script_dir+'/variantPipeline.bpipe.groovy',output_dir)

# add variables to the bpipe config file to pass them to the pipeline
with open(output_dir+'/variantPipeline.bpipe.config.groovy','w') as config:
    config.write('REFERENCE='+'\"'+ ref+ '\"'+'\n') # The name of the reference files for bowtie alignment wit
    config.write('REFERENCE_FA='+ '\"'+ref+ '.fa' '\"'+'\n') # The reference file fasta to be used in the deepSNV step relative to segment and locations to call variants
    config.write('SCRIPTS='+ '\"'+script_dir+ '\"'+'\n') # The scripts dir 
    config.write('LIBRARY_LOCATION='+ '\"'+lib_dir+'\"'+ '\n') # The library dir
    config.write('CONTROL='+ '\"'+control+ '\"'+'\n') # The name of the plasmid control
    config.write('DISP='+ '\"'+disp+ '\"'+'\n')# The Dispersion estimation to be used
    config.write('P_CUT='+ '\"'+p_cut+ '\"'+'\n') # The p cut off
    config.write('P_COM_METH='+ '\"'+method+ '\"'+'\n') # The combination method used to combine the pvalues from each strand
    config.write('INPUT_DIR='+ '\"'+input_dir+ '\"'+'\n') # copy the input dir to the config file to help find the control when running in bam
    config.write('OR='+ '\"'+open_reading+ '\"'+'\n') # copy the open reading frame file
#throttled to 3 processors to be a good neighbor.
#note that running unthrottled can result in errors when bpipe overallocates threads/memory

# if bam==True: #If bam is set only look for the bam files
#     if test==False:
#         command= bpipe_command + " run -r " + output_dir +  "/variantPipeline.postalign.bpipe.groovy " + input_dir + "/*.bam"
#     else:
#         command=bpipe_command + " test " + output_dir +  "/variantPipeline.postalign.bpipe.groovy " + input_dir +"/*.bam"
# else: # Otherwise start with the fastqs
if test==False:
	command= bpipe_command + " run -r " + output_dir +  "/variantPipeline.bpipe.groovy " + input_dir + "/*.fastq"
else:
	command=bpipe_command + " test " + output_dir +  "/variantPipeline.bpipe.groovy " + input_dir +"/*.fastq"
print "submitting command: \n"+command



s.call(command,shell=True) # rub bpipe command

sys.exit(0)
