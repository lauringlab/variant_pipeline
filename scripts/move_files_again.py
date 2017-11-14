import argparse
import shutil
import os
import errno
import re
parser = argparse.ArgumentParser(description='So you\'ve run your analysis and deleted the copied fastq files to save space and now you need to run it again. This script takes in the renaming log file in the form initialfile\ttext\tfinal file and recopies the files. \n This is also useful in prepping to upload files to the SRA. by activating -sra_dir and -sra_file you can specify the samples and final directory for files you need to upload. This option appends the original run directory to the sample name for clarity.')
parser.add_argument("s",nargs='+',metavar='s',help='The renaming file.')
parser.add_argument("-r",action="store_true",dest='run', default=False,help = "Boolean switch to actually run the program.")
parser.add_argument("-sra",action="store_true",dest='sra', default=False,help = "Boolean switch to activate SRA mode.")
parser.add_argument("-sra_dir",action="store",dest='sra_dir',help = "The directory that will hold the SRA ready files")
parser.add_argument("-sra_files",action="store",dest='sra_files',help = "A text file containing the sample IDs, (one per line). Final file names that include the sample ID will be copied")
parser.add_argument("-run_name",action="store",dest='run_name',help = "Run name to add to sra file")

args=parser.parse_args()
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def recopy(rename_log,run):
    with open(rename_log,'r') as file:
        for line in file:
            line=line.split("\t")
            old=line[0].strip()
            new=line[2].strip()
            #print new 
            if run==True:
                shutil.copy2(old,new)
            else:
                print "test mode"
            print "copying %s to %s" % (old,new) 

def get_samples(sample_txt):
    samples = []
    with open(sample_txt,'r') as file:
        for line in file:
            line=line.strip()
            samples.append(line)
    return(samples)

def sracopy(rename_log,run,sra_dir,sra_files,run_name):
    base_name = re.compile(r'(.*)(\.fastq.*$)')
    samples = get_samples(sra_files)
    with open(rename_log,'r') as file:
        for line in file:
            line=line.split("\t")
            old=line[0].strip()
            new=line[2].strip()
            # adjust new to go to sra directory    
            file_name = new.split("/")[-1]
            #dir_name = new.split("/")[-2]
            
            final_new = "%s/%s" %(sra_dir, base_name.sub("\\1.%s\\2" %run_name,file_name))
            #final_new = sra_dir+"/"re.sub('\.fastq.*$', , somestring)+file_name+"_"+dir_name

            for sample in samples:
                    if sample in new:
                        print "sample: %s, new %s" %(sample,new)
                        if run==True:
                            shutil.copy2(old,final_new)
                        else:
                            print "test mode"
                        print "copying %s to %s" % (old,final_new)

def main():
    if args.sra : 
        make_sure_path_exists(args.sra_dir) # Make sra directory if needed
        sracopy(args.s[0],args.run,args.sra_dir,args.sra_files,args.run_name)	

    else:
        recopy(args.s[0],args.run)
main()
