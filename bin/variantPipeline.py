import os
import os.path
import argparse
import shutil
import subprocess as s
import sys



parser = argparse.ArgumentParser(description='This is a wrapper to set up and run the bpipe command')

parser.add_argument('-i',action='store',dest='input_dir',help='Directory containing the input fastqs')
parser.add_argument('-o',action='store',dest='output_dir',help='The final directory that will hold the output. If it does\'t exsit it will be made')
parser.add_argument('-r',action='store',dest='ref',help='The name of the reference files used for bowtie alignment')
parser.add_argument('-p',action='store',dest='control',help='The sample name of the plasmid control used for varinat calling')
parser.add_argument('-t',action='store_true',dest='test',default=False,help='Boolean switch to run program in test mode. Everything will be set up but bpipe will run in test mode')

parser.add_argument('-d',action='store',dest='disp',help='Dispersion estimation to be used in deepSNV. options are c("two.sided","one.sided","bin"')

args=parser.parse_args()

input_dir=os.path.abspath(args.input_dir)
output_dir=os.path.abspath(args.output_dir)
ref=os.path.abspath(args.ref)
control=args.control
disp=args.disp
bin_dir=os.path.dirname(os.path.realpath(__file__))
script_dir=os.path.abspath(bin_dir+'/..'+'/scripts/')
lib_dir=os.path.abspath(bin_dir+'/..'+'/lib/')
bpipe_command=lib_dir+'/bpipe-0.9.8.7/bin/bpipe'
test=args.test

print "Processing fastqs from " + input_dir
print "Results will be saved to " + output_dir
print "Using " + ref +" for a reference and \n" + control + " as the control sample"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.chdir(output_dir)

shutil.copy(script_dir+'/variantPipeline.bpipe.stages.groovy',output_dir)
shutil.copy(script_dir+'/variantPipeline.bpipe.groovy',output_dir)


# add variables to config reference to config file
with open(output_dir+'/variantPipeline.bpipe.config.groovy','w') as config:
    config.write('REFERENCE='+'\"'+ ref+ '\"'+'\n')
    config.write('REFERENCE_FA='+ '\"'+ref+ '.fa' '\"'+'\n')
    config.write('SCRIPTS='+ '\"'+script_dir+ '\"'+'\n')
    config.write('LIBRARY_LOCATION='+ '\"'+lib_dir+'\"'+ '\n')
    config.write('CONTROL='+ '\"'+control+ '\"'+'\n')
    config.write('DISP='+ '\"'+disp+ '\"'+'\n')

#throttled to 8 processors to be a good neighbor.
#note that running unthrottled can result in errors when bpipe overallocates threads/memory
if test==False:
    command= bpipe_command + " run -n 8 -r " + output_dir +  "/variantPipeline.bpipe.groovy " + input_dir + "/*.fastq"
else:
    command=bpipe_command + " test -n 8 " + output_dir +  "/variantPipeline.bpipe.groovy " + input_dir +"/*.fastq"
print "submitting command: \n"+command

#p = s.Popen("time", stdin=s.PIPE, stdout=s.PIPE, stderr=s.STDOUT)
#p.stdin.write(command)
#out, err = p.communicate()

s.call(command,shell=True)#, stdout=subprocess.PIPE)
#output = process.communicate()[0]

sys.exit(0)
