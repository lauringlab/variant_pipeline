import os
import os.path
import argparse
import shutil
import subprocess as s
import yaml
import sys

parser = argparse.ArgumentParser(description='This is a wrapper to set up and run the bpipe command')


parser.add_argument('pipeline',metavar = 'pipeline',nargs='+',
					help = 'the name of the bpipe pipeline script to be run')
parser.add_argument('input_files',metavar = 'input_files',nargs='+',
					help = 'the relative path and wildcard to input files')
parser.add_argument('output_dir',metavar = 'output_dir',nargs='+',
					help = 'the directory to hold the outputs')
parser.add_argument('options', metavar='options', nargs='+',
                    help='A YAML options file mimicing the one found in the bin directions')
                    
parser.add_argument('-t',action='store_true',dest='test',default=False,
					help='Boolean switch to run program in test mode. Everything will be set up but bpipe will run in test mode')

args=parser.parse_args()

with open(args.options[0], 'r') as stream:
    try:
        options=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise("YAML error")


## Give the input arguments better names ##

options_file = os.path.abspath(args.options[0])
pipeline = os.path.abspath(args.pipeline[0])

# get absolute path to input files
input_file_string = args.input_files[0]
input_path = ""
input_file_split = input_file_string.split("/")
i = 0
while i <(len(input_file_split)-1):
    input_path = input_path+input_file_split[i]+"/"
    i+=1

input_files_abspath=os.path.abspath(input_path)
input_files = str(input_files_abspath)+"/"+input_file_split[-1]
output_dir=os.path.abspath(args.output_dir[0])


# options that need absolute paths
path_options = ["REFERENCE","REFERENCE_FA","OPEN_READING"]

# 
#     
# ref=os.path.abspath(options["ref"])
# control=options["control"]
# disp=options["disp"]
# p_cut=options["p"]
# method=options["method"]
# 
# ## options for processing variants ###
# open_reading=os.path.abspath(options["open_reading"])
# stringent_freq = options["stringent_freq"]

bin_dir=os.path.dirname(os.path.realpath(__file__)) # The path to this file so we can find the scripts and lib
script_dir=os.path.abspath(bin_dir+'/..'+'/scripts/')# The path to the scripts dir relative to this location
lib_dir=os.path.abspath(bin_dir+'/..'+'/lib/') # The path to the lib dir relative to this location
bpipe_command=lib_dir+'/bpipe-0.9.8.7/bin/bpipe' # The path to the bpipe command relative the lib dir.
test=args.test
#bam=args.bam

print("Processing fastqs from " + input_files)
print("Results will be saved to " + output_dir)
#print("Using " + ref +" for a reference and \n" + control + " as the control sample")


cwd = os.getcwd()
os.chdir(bin_dir+"/..")
git_command="git rev-list HEAD |head -n 1"
version=s.check_output(git_command,shell=True)
os.chdir(cwd)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

local_pipeline_name = os.path.basename(pipeline)
config_file = output_dir+'/'+local_pipeline_name+'.config.groovy'
local_pipeline = output_dir+'/'+local_pipeline_name
print local_pipeline
with open(local_pipeline,'w') as local_p:
    local_p.write("// Loading config file \n")
    local_p.write("load \'"+config_file+"\' \n")
    with open(pipeline,"r") as original_p:
        for line in original_p:
            local_p.write(line)



# add variables to the bpipe config file to pass them to the pipeline
with open(config_file,'w') as config:
    config.write('//This pipeline was run on with commit : '+ version.decode() +'\n')
    config.write('LIBRARY_LOCATION='+ '\"'+lib_dir+'\"'+ '\n') # The library dir
    config.write('SCRIPTS='+ '\"'+script_dir+ '\"'+'\n') # The scripts dir 

    for option in options:
        opt = options[option]
        if option in path_options:
            opt = os.path.abspath(opt)
        config.write(option+'='+'\"' + str(opt) + '\"'+'\n')




## If the output dir does not exist make it 

os.chdir(output_dir)

# Copy the stages script and the pipeline to the output dir
shutil.copy(script_dir+'/variantPipeline.bpipe.stages.groovy',output_dir)

# copy pipeline script to local file adding the load command
#shutil.copy(pipeline,output_dir)


        
#     config.write('REFERENCE='+'\"'+ ref+ '\"'+'\n') # The name of the reference files for bowtie alignment wit
#     config.write('REFERENCE_FA='+ '\"'+ref+ '.fa' '\"'+'\n') # The reference file fasta to be used in the deepSNV step relative to segment and locations to call variants
#     config.write('SCRIPTS='+ '\"'+script_dir+ '\"'+'\n') # The scripts dir 
#     config.write('CONTROL='+ '\"'+control+ '\"'+'\n') # The name of the plasmid control
#     config.write('DISP='+ '\"'+disp+ '\"'+'\n')# The Dispersion estimation to be used
#     config.write('P_CUT='+ '\"'+str(p_cut)+ '\"'+'\n') # The p cut off
#     config.write('P_COM_METH='+ '\"'+method+ '\"'+'\n') # The combination method used to combine the pvalues from each strand
#     config.write('OR='+ '\"'+open_reading+ '\"'+'\n') # copy the open reading frame file
#     config.write('OPTIONS=' + '\"' + options_file + '\"\n') # Copy the options file 
#     config.write('STRINGENT_FREQ=' + '\"' + str(stringent_freq) + '\"\n') # The frequency below which deepSNV and stringent thresholds are needed
#note that running unthrottled can result in errors when bpipe overallocates threads/memory

if test==False:
    command= bpipe_command + " run -r " +  local_pipeline + " " + input_files
else:
	command=bpipe_command + " test " + local_pipeline + " " + input_files
print("submitting command: \n"+command)



s.call(command,shell=True) # rub bpipe command

sys.exit(0)
