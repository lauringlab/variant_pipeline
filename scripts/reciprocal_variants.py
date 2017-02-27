from __future__ import division
import argparse
import pandas as pd
import pysam
import copy
import numpy as np
import yaml 


parser = argparse.ArgumentParser(description='This script is designed to read in a csv of putative variant calls from deepSNV and query the bam file for the number of reads matching the reference base at that position.',usage ="python reciprocal_variants.py sample.bam sample.csv")

parser.add_argument('bam', metavar='bam', nargs='+',
                    help='The bam file of the sample')
parser.add_argument('csv', metavar='csv', nargs='+',
                    help='The deepSNV csv file')
parser.add_argument('out_csv', metavar='out_csv', nargs='+',
                    help='The csv file that will contain the output')
parser.add_argument('options', metavar='options', nargs='+',
                    help='A YAML options file mimicing the one found in the bin directions')


def var_setup(row):
    df=copy.deepcopy(row)
    df["mapq"]=[]
    df["phred"]=[]
    df["read_pos"]=[]
    return df
def ref_setup(row):
    df=var_setup(row)
    df["var"]=df["ref"]
    nonapplicable_items=["p.val","sigma2.freq.var","n.ctrl.fw","n.ctrl.bw","raw.p.val"]
    zero_items=["n.tst.fw","n.tst.bw","freq.var"]
    df[nonapplicable_items]="NA"
    df[zero_items]=0
    df["mutation"]=df["chr"]+"_"+df["ref"]+str(df["pos"])+df["var"]
    return df

def update_base(read,df,n_counts):
    if  not read.is_del and not read.is_refskip: 
        called_base=read.alignment.query_sequence[read.query_position] # what is the called base at the position
        called_phred=read.alignment.query_qualities[read.query_position]
        if called_phred>30 and called_base==df["var"]:
            df["mapq"].append(read.alignment.mapping_quality)
            df["phred"].append(called_phred)
            df["read_pos"].append(read.query_position)
            if n_counts:
                if read.alignment.is_reverse:
                    df["n.tst.bw"]+=1
                elif not read.alignment.is_reverse:
                    df["n.tst.fw"]+=1
    return df

def finalize(df):
    df["MapQ"]=np.mean(df.mapq)
    df["Phred"]=np.mean(df.phred)
    df["Read_pos"]=np.mean(df.read_pos)
    df=df.drop(["mapq","phred","read_pos"])
    return df


args=parser.parse_args()    


with open(args.options[0], 'r') as stream:
    try:
        options=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise "YAML error"


bam= pysam.AlignmentFile(args.bam[0], "rb") 
variants=pd.read_csv(args.csv[0],index_col=False)
out_csv=args.out_csv[0]

def get_reference(row,bam):
    # grab some data from the first series.
    var=var_setup(row)
    ref=ref_setup(row)
    py_pos=var.pos-1 # pythonify the position
    for pileupcolumn in bam.pileup(var.chr,py_pos,py_pos+1,truncate=True,stepper="all",max_depth=1E6): # look at the range py_pos to py_pos+1 and all the bases that align to those positions. Think of a column on each position like the samtools output of mpileup
        if not pileupcolumn.pos==py_pos: 
            raise RuntimeError('This is not the correct position')
        for pileupread in pileupcolumn.pileups: #Now for each read in that column
            var=update_base(pileupread,var,n_counts=False)
            ref=update_base(pileupread,ref,n_counts=True)
    # verify everything is up to snuff
    
    if not len(var.mapq)==(var["n.tst.fw"]+var["n.tst.bw"]):
        raise RuntimeError('We did not find the same number of variant bases as deepSNV - something is up')
    if not len(ref.mapq)==(ref["n.tst.fw"]+ref["n.tst.bw"]):
        raise RuntimeError("Didn't get the same number of reference qualities as bases- something is up")

    # Finalize series
    var=finalize(var)
    ref=finalize(ref)
    ref["freq.var"]=(ref["n.tst.fw"]+ref["n.tst.bw"])/(ref["cov.tst.fw"]+ref["cov.tst.bw"])
    return([var,ref])
    
inferred_qual=pd.DataFrame()
for index, row in variants.iterrows():
    meta_bases=get_reference(row,bam)                                                    
    inferred_qual=inferred_qual.append(meta_bases[0])
    inferred_qual=inferred_qual.append(meta_bases[1])

inferred_qual_nodup=inferred_qual.drop_duplicates() # If there are 2 variants at any position we will grab the reference data twice. This removes those duplicate rows.

if not options["infer"]: # infer is a boolean - False if we don't want to infer the variants.
    no_infer=inferred_qual_nodup.loc[ (inferred_qual_nodup['ref']!=inferred_qual_nodup['var'])]
    no_infer.to_csv(out_csv)
else:
    inferred_qual_nodup.to_csv(out_csv)



