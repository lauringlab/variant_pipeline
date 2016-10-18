#!/Users/jt/.virtualenvs/sci-py2.7/bin/python

import os
import os.path
import argparse
import shutil
import subprocess as s
import sys



parser = argparse.ArgumentParser(description='This is a wrapper to set up and run the bpipe pipeline that processes the deepSNV fasta files')

parser.add_argument(action='store',dest='input_dir',nargs='+',help='Directory containing the input fasta')
#parser.add_argument(action='store',dest='coding',nargs='+',help='The name of the reference fasta file for aligning the sequences')
parser.add_argument(action='store',dest='ref',nargs='+',help='The name of the reference fasta file that is a template for the parsing')
parser.add_argument('-t',action='store_true',dest='test',default=False,help='Boolean switch to run program in test mode. Everything will be set up but bpipe will run in test mode')

args=parser.parse_args()

## Give the input arguments better names ##
input_dir=os.path.abspath(args.input_dir[0])
main_data_dir=input_dir+"/.."
ref=os.path.abspath(args.ref[0])
#coding=os.path.abspath(args.coding[0])
script_dir=os.path.dirname(os.path.realpath(__file__)) # The path to this file so we can find the scripts and lib
bpipe_command='~/variant_pipeline/lib/bpipe-0.9.8.7/bin/bpipe' # The path to the bpipe command relative the lib dir.
test=args.test


print "Processing fastas from " + input_dir


# add variables to the bpipe config file to pass them to the pipeline
with open('./fasta.config.groovy','w') as config:
    config.write('REFERENCE_FA='+'\"'+ ref+ '\"'+'\n') # The name of the reference for deconcatination
#    config.write('CODING_FA='+'\"'+ coding+ '\"'+'\n') # The name of the reference csv for deconcatination
    config.write('SCRIPTS='+ '\"'+script_dir+ '\"'+'\n') # The scripts dir
    config.write('MAIN_DIR='+ '\"'+main_data_dir+ '\"'+'\n') # copy the input dir to the config file to help find the control when running in bam
#throttled to 4 processors to be a good neighbor.
#note that running unthrottled can result in errors when bpipe overallocates threads/memory


if test==False:
    command= bpipe_command + " run -n 4 -r " + script_dir+ "/consensus.pipe.groovy " + input_dir + "/*.fa"
else:
    command = bpipe_command + " test -n 4 " + script_dir +  "/consensus.pipe.groovy " + input_dir +"/*.fa"
print("submitting command: \n"+command)
s.call(command,shell=True) # rub bpipe command

sys.exit(0)
