from __future__ import division
import pysam
import numpy as np
import yaml 


parser = argparse.ArgumentParser(description='This scipts takes a bam file and identifies the consensus file sequence of the sample.',usage ="python consensus.py sample.bam sample.fa")

parser.add_argument('bam', metavar='bam', nargs='+',
                    help='The bam file of the sample')
parser.add_argument('sample.fa', metavar='fasta', nargs='+',
                    help='the consensus sequence output file')

class locus(object):
    """ A base which has the following characteristics 
        Atrributes :
            Chr : chromosome
            pos : position on chromosome
            counts : a dictionary of { A : the number of As
            T : the number of Ts
            C : the number of Cs
            G : the number of Gs}
        methods : 
            update : update counts
            calc_freqs :  calculate the frequency of each base
            consensus : caculate the consensus sequence as this position
    """
    def __init__(self,chr,pos):
        """return a base object with chr and pos and starting counts of 0
        Posisition is base 1 like most biologists are use to."""
        self.chr = chr
        self.pos = pos
        self.counts = {'A':0,'T':0,'C':0,'G':0,'-':0}
        self.freqs = {'A':0,'T':0,'C':0,'G':0,'-':0}
        self.coverage = 0
    def update(self,base):
        self.coverage = self.coverage+1
        self.counts[base] = self.counts[base]+1
    def calc_freqs(self):
        for base in self.counts.keys():
            self.freqs[base] = self.counts[base]/self.coverage
    def consensus(self,cutoff):
        self.calc_freqs()
        consensus = {k:v for (k,v) in self.freqs.items() if v > cutoff}
        if len(consensus)==1:
            return list(consensus.keys())[0]
        else:
            return "N"
        

class segment(object):
    """ A sequence like object made up of base objects
            Attriutes:
                chr - the name of the chr
                seq - a list of loci
            methods:
                append: append another position
                update: update base count at a position (loci base 1)
                consensus: calcuate the consensus sequence"""
    def __init__(self,chr):
        self.chr = chr
        self.seq = []
    def append_loci(self,loci):
        if type(loci) is not locus:
            raise ValueError('Only class locus can be appended to a segment object')
        if loci.chr!=self.chr:
            raise ValueError('The loci chr does not match the segement chr')
        if loci.pos!=(len(self.seq)+1):
            raise ValueError('The position of the loci does not match current segment length')
        self.seq.append(loci)
    def consensus(self,cutoff):
        seg_consensus = ""
        for loci in self.seq:
            seg_consensus=seg_consensus+loci.consensus(cutoff)
        return seg_consensus
    def locus(self,pos):
        loci = [x for x in self.seq if x.pos==pos]
        return loci[0]


class genome(object):
            

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
    if len(df["mapq"])==0:
        df["MapQ"]="NA"
        df["Phred"]="NA"
        df["Read_pos"]="NA"
    else:
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
        raise("YAML error")


bam= pysam.AlignmentFile(args.bam[0], "rb") 
variants=pd.read_csv(args.csv[0],index_col=False)
out_csv=args.out_csv[0]

def get_reference(row,bam,condition):
    # grab some data from the first series.
    if condition == "var":
        df=var_setup(row)
        nc = False # see update base command - we don't add to the counts here
    elif condition =="ref":
        df=ref_setup(row)
        nc = True # see update base command - we do add to the counts here
    else:
        raise RuntimeError('Please enter ref or var as the condition.')
    
    py_pos=df.pos-1 # pythonify the position
    for pileupcolumn in bam.pileup(df.chr,py_pos,py_pos+1,truncate=True,stepper="all",max_depth=1E6): # look at the range py_pos to py_pos+1 and all the bases that align to those positions. Think of a column on each position like the samtools output of mpileup
        if not pileupcolumn.pos==py_pos: 
            raise RuntimeError('This is not the correct position')
        for pileupread in pileupcolumn.pileups: #Now for each read in that column
            df=update_base(pileupread,df,n_counts=nc)
    # verify everything is up to snuff
    
    if condition=="var" and  len(df.mapq)!=(df["n.tst.fw"]+df["n.tst.bw"]):
        print(df)
        raise RuntimeError('We did not find the same number of variant bases as deepSNV - something is up')
    if not len(df.mapq)==(df["n.tst.fw"]+df["n.tst.bw"]):
        print(df)
        raise RuntimeError("Didn't get the same number of reference qualities as bases- something is up")

    # Finalize series
    df=finalize(df)
    if condition == "ref":
        df["freq.var"]=(df["n.tst.fw"]+df["n.tst.bw"])/(df["cov.tst.fw"]+df["cov.tst.bw"])
    return(df)
    
inferred_qual=pd.DataFrame()
no_inferred_qual=pd.DataFrame()

for index, row in variants.iterrows():
    
    if row["freq.var"]>=options["stringent_freq"] and row["ref"]==row["var"]: # if the allele is not subject to stringent measures and it is the reference base we move on saving it for later.
        no_inferred_qual =  no_inferred_qual.append(row,ignore_index=True)
    elif row["freq.var"]>=options["stringent_freq"] and row["ref"]!=row["var"]: # if the allele is not subject to stringent measures and it differs from the reference base then we look for the reference base in the stringent fraction as deepSNV.R would not have looked there 
        no_inferred_qual =  no_inferred_qual.append(row,ignore_index=True)
        ref=get_reference(row,bam,"ref")
        if ref["freq.var"]<options["stringent_freq"] and ref["freq.var"]>0: # ensures we only add it if the frequency puts it in the stringent fractor.
           inferred_qual =  inferred_qual.append(ref,ignore_index=True)
    elif row["freq.var"]<options["stringent_freq"]:
        var = get_reference(row,bam,"var")
        inferred_qual=inferred_qual.append(var,ignore_index=True)
### Add meta data columns to rows that don't need them.
print(no_inferred_qual.head())
print(inferred_qual.head())

no_inferred_qual["MapQ"]='NA'
no_inferred_qual["Phred"]='NA'
no_inferred_qual["Read_pos"]='NA'

## Combine rows with meta data and rows without it.

all_sites = inferred_qual.append(no_inferred_qual,ignore_index=True)

if not options["infer"]: # infer is a boolean - False if we don't want to infer the variants.
    no_infer=all_sites.loc[ (all_sites['ref']!=all_sites['var'])]
    no_infer.to_csv(out_csv)
else:
    all_sites.to_csv(out_csv)



