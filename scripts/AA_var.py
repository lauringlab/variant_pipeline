from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import pandas as pd
import os
import tempfile
import sys
import subprocess
import numpy
import copy
import argparse
import yaml
from fasta_functions import *

#    in_csv = ["~/Documents/Analysis/scratch/Filter_var/2377_1.removed.mapq.sum.filtered.csv"]
#    OR = ["/Users/jt/Documents/Analysis/scratch/polio_multiple_OR.fa"]
#    Ref = ["/Users/jt/Documents/Analysis/scratch/polio_multiple.fa"]
parser = argparse.ArgumentParser(description='This script takes a consensus fasta, a fasta file of open reading frames (whose names include the segment names in the consensus fasta), and a csv of variant calls as outputed from deepSNV. It then identifies mutations as nonsynonymous or synonymous. It currently requires that the open reading frame and the reference sequence have the same number of amino acids. Any difference will through an error.')

parser.add_argument('Ref', metavar='Ref', nargs='+',
                    help='The reference consensus fasta file. Variants will be added to this background')

parser.add_argument('in_csv', metavar='in_csv', nargs='+',
                    help='The variant csv.')
parser.add_argument('out_csv', metavar='out_csv', nargs='+',
                    help='The output csv.')
parser.add_argument('options', metavar='options', nargs='+',
                    help='A YAML options file mimicing the one found in the bin directions')

opts = parser.parse_args()

# I am making the options a class since that is how the code was structured prior to reading in options from a Yaml. The yaml could include the class in the future
with open(opts.options[0], 'r') as stream:
    try:
        options=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise "YAML error"


class test_args:
    Ref = [opts.Ref[0]]
    in_csv= [opts.in_csv[0]]
    out_csv= [opts.out_csv[0]]
    muscle_path = options["muscle_path"]
    OR = options["open_reading"]
    classification = options["classification"]
    
args = test_args()
    
def no_nones(x):
    return [y for y in x if y is not None]

def trim_noncoding(series):
    # If the variant is never in an open reading frame return [noncoding] for coding_pos Ref_AA AA_pos and 
    
    if sum([x is None for x in series["coding_pos"]]) == len(series["coding_pos"]): # is true if all entries are None
        series["coding_pos"]=['Noncoding']
        series["Ref_AA"]=['Noncoding']
        series["AA_pos"]=['Noncoding']
        series["Class"]=['Noncoding']
        series["Var_AA"]=["Noncoding"]
        series["OR"]=["Noncoding"]
    else: # remove the None's from the lists for variants lie in at least one open reading frame
        series["OR"]=no_nones(series["OR"])
        series["coding_pos"]=no_nones(series["coding_pos"])
        series["Ref_AA"]=no_nones(series["Ref_AA"])
        series["AA_pos"]=no_nones(series["AA_pos"])
        series["Class"]=no_nones(series["Class"])
        series["Var_AA"]=no_nones(series["Var_AA"])
    return(series)

    


#class test_args:
#    in_csv = ["~/Documents/Analysis/scratch/Filter_var/2377_1.removed.mapq.sum.filtered.csv"]
#    OR = ["/Users/jt/Documents/Analysis/scratch/polio_multiple_OR.fa"]
#    Ref = ["/Users/jt/Documents/Analysis/scratch/polio_multiple.fa"]
#args=test_args()
variants=pd.read_csv(args.in_csv[0],index_col=0)

# setting up empty lists - OR - Open reading frame
###################################################################
#
#
#
#             The good stuff                                         
#
#
#
#
####################################################################


variants["OR"]= [list() for x in range(len(variants.index))]
variants["coding_pos"]= [list() for x in range(len(variants.index))]
variants["Ref_AA"]= [list() for x in range(len(variants.index))]
variants["AA_pos"]= [list() for x in range(len(variants.index))]
variants["Var_AA"]= [list() for x in range(len(variants.index))]
variants["Class"]= [list() for x in range(len(variants.index))]



coding=ReadFASTA(args.OR)
reference = ReadFASTA(args.Ref[0])

for ref in reference:


    for code in coding:
        seg = ref.id
        OR= code.id
        if seg in OR:
            print "Working with variants on %s" % seg
            seg_var = variants.index[variants["chr"]==seg]
            fixed=variants.loc[(variants["freq.var"]>0.5) & (variants["chr"]==seg)] # anything above 50% will be in the consensus or "fixed"
            # Mutate the consensus variants to make the plasmid consensus into the sample consensus.
            fixed_ref=mutate(ref,fixed)
            # Align
            fix_ref_coding=Align([fixed_ref,code],args.muscle_path)
            ref_coding=Align([ref,code],args.muscle_path)
            ref_coding_deep=copy.deepcopy(ref_coding)
            # Trim to just coding sequence
            ref_trimmed=StripGapsToFirstSequence([ref_coding[1],ref_coding[0]])
            fixed_ref_trimmed = StripGapsToFirstSequence([fix_ref_coding[1],fix_ref_coding[0]])
            # Translate to amino acid sequence
            # remove any gaps in ref seq 
            if fixed_ref_trimmed.seq.count("-")>0:
                raise ValueError, "Gaps were found in the open reading frame - check sequence for insertions" # think about writing in code to tranlate these if they are inframe
            if args.classification == 'control':
                ref_trans=ref_trimmed.seq.translate() # This is reference that is used to call the reference AA
            elif args.classification=='sample':
                ref_trans = fixed_ref_trimmed.seq.translate()
            else:
                raise "unknown classification reference - please use either 'sample' or 'control' only."
            # Set fixed mutations at the nucleotide level
            # Get coding pos for appropriate variants
            # Note GetCorrespondingResidue and the AA_pos are in base 1 not the usual base 0 python -  I want to fix this at some point
            for i in seg_var:
                variants.loc[i,"coding_pos"].append(GetCorrespondingResidue(ref_coding_deep,variants.loc[i,"pos"]))
            # Get AA_pos for appropirate variants
                variants.loc[i,"AA_pos"].append( catch(lambda : (variants.loc[i,"coding_pos"][-1]-1)//3+1 )) # the [-1] is to get the most recently added one
            # Get reference AA ""
                variants.loc[i,"Ref_AA"].append(catch(lambda : ref_trans[int(variants.loc[i,"AA_pos"][-1])-1])) # the second minus 1 is to convert to python numbering
            
            # Cycle through each variants
                mutation = var_aa(fixed_ref_trimmed,variants,i)
                variants.loc[i,"Var_AA"].append(mutation["AA"])
                variants.loc[i,"Class"].append(mutation["Class"])
            # Add open reading frame
                variants.loc[i,"OR"].append(OR)

variants=variants.apply(trim_noncoding,axis=1)


variants.to_csv(args.out_csv[0])

    
    
    
    
    
    

    
    

