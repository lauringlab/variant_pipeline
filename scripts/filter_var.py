import pandas as pd
import numpy as np
import argparse
import yaml

parser = argparse.ArgumentParser(description='The purpose of this script is to take in variant calls from deepSNV and filter them based on a set of quality scores (non of them required. Also it can add metadata if a metadata file is provided. The script can also infer minor variants for positions that differ at the consensus level from the plasmid control.',usage ="python scripts/filter_var.py data/process/HK_1/Variants/all.sum.csv ./test.csv -mapping 30 -phred 35 -p 0.01 -freq 0.01 -pos 32 94  -run HK_6 -meta ./data/reference/all_meta.csv -infer")


parser.add_argument('variants',metavar='variants',nargs='+',
			help = 'the input csv file')
parser.add_argument('out_csv',metavar='out_csv',nargs='+',
                        help = 'the out_put csv file')
parser.add_argument('options', metavar='options', nargs='+',
                    help='A YAML options file mimicing the one found in the bin directions')
opts = parser.parse_args()

# I am making the options a class since that is how the code was structured prior to reading in options from a Yaml. The yaml could include the class in the future
with open(opts.options[0], 'r') as stream:
    try:
        options=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise "YAML error"


variants = opts.variants[0]
phred = options["phred"]
mq = options["mapping"]
freq = options['freq']
pos = options['pos']
pval = options['p_cut']
run= options['run']
stringent_freq = options['stringent_freq']
out_csv= opts.out_csv[0]



def filter(x,mq,phred,freq,stringent_freq,pval,pos):
    stringent_variants = x.loc[ (x['ref']!=x['var']) & (x['MapQ']>mq) & (x['Phred']>phred) & (x['freq.var']>freq) & (x['freq.var']<stringent_freq) & (x['p.val']<pval ) & (x['Read_pos']>pos[0]) & (x['Read_pos']<pos[1])]
    not_stringent = x.loc[(x['freq.var']>stringent_freq)]
    stringent_ref  = x.loc[ (x['ref']==x['var']) & (x['MapQ']>mq) & (x['Phred']>phred) & (x['freq.var']>freq) & (x['freq.var']<stringent_freq) & (x['Read_pos']>pos[0]) & (x['Read_pos']<pos[1])] # same as above but no p.val
    out = stringent_variants.append(not_stringent,ignore_index=True)
    out=out.append(stringent_ref,ignore_index=True)
    out=out.drop_duplicates()
    return(out)
#print mq

# read in the csv and apply the subsetting function
#print pos
data=pd.DataFrame.from_csv(variants,index_col=None)
data.Id.apply(str)
out = filter(data,mq,phred,freq,stringent_freq,pval,pos)

if out.shape[0]>0: # if there are some variants left
    if run != None:
        print("Adding run label")
        out.loc[:,'run']=run
else:
    print("There were no variants left - There may have been none to begin with")
    
print("Saving output as " + out_csv)
out.to_csv(out_csv)







