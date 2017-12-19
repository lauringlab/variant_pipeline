import os
import os.path
import argparse
import shutil
import subprocess as s
import sys
import pandas as pd


parser = argparse.ArgumentParser(description='This is a wrapper to set up and run the bpipe pipeline that downloads fastq file from the SRA and names them accroding to a csv file that contains at least the headings sra_accession","sample","fastq1","fastq2" were fastq* contains the path and name of the final fastq files.')

parser.add_argument(action='store',dest='csv_in',nargs='+',help=' a csv containing the files to be downloaded and moved.')
parser.add_argument(action='store',dest='branches',nargs='+',help='The name of the brances file that will used to run the download pipe')
parser.add_argument(action='store',dest='out_dir',nargs='+',help='The directory that will hold the final outputs')
parser.add_argument('-p',action='store',dest = 'processors',default=4,type=int,help = 'the number of processors to use - default is 4')
parser.add_argument('-t',action='store_true',dest='test',default=False,help='Boolean switch to run program in test mode. Everything will be set up but bpipe will run in test mode')

args=parser.parse_args()

# Get the paths to key files relative to this one

bin_dir=os.path.dirname(os.path.realpath(__file__)) # The path to this file so we can find the scripts and lib
script_dir=os.path.abspath(bin_dir+'/..'+'/scripts/')# The path to the scripts dir relative to this location
lib_dir=os.path.abspath(bin_dir+'/..'+'/lib/') # The path to the lib dir relative to this location
bpipe_command=lib_dir+'/bpipe-0.9.8.7/bin/bpipe' # The path to the bpipe command relative the lib dir.
test=args.test
out_dir = os.path.abspath(args.out_dir[0])

print("downloading files from %s" % args.csv_in[0])

# Check that ouput directories exist

csv_file = pd.read_csv(args.csv_in[0])

if not os.path.exists(out_dir):
    print("Directory %s does not exist. making it now" %d)
    os.makedirs(d)

# make branches file 
# This is useful in case the input csv has some columns we don't need
print(args.branches[0])
# Copy the stages script and the pipeline to the output dir
with open(args.branches[0],'w') as b:
    for index,row in csv_file.iterrows():
        b.write("%s\t%s\t%s\t%s\n" %(row['sample'],row['sra_accession'],row['fastq1'],row['fastq2']))

# add variables to the bpipe config file to pass them to the pipeline
with open('./download.config.groovy','w') as config:
    config.write('BRANCHES ='+ '\"'+os.path.abspath(args.branches[0])+ '\"'+'\n')
    config.write('OUT ='+ '\"'+ out_dir + '\"'+'\n')
# copy bpipe script to the working directory
shutil.copy(script_dir+'/download.fastq.groovy',"./")

#throttled to 4 processors to be a good neighbor.
#note that running unthrottled can result in errors when bpipe overallocates threads/memory


if test==False:
    command= "%s run -n %d -r ./download.fastq.groovy" %(bpipe_command, args.processors)
else:
    command = "%s test -n %d -r ./download.fastq.groovy" %(bpipe_command, args.processors)
print("submitting command: \n"+command)
s.call(command,shell=True) # rub bpipe command

sys.exit(0)

