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


class test_args:
    variants = [opts.variants[0]]
    phred = options["phred"]
    mq = options["mapping"]
    freq = options['freq']
    pos = options['pos']
    pval = options['p_cut']
    run= options['run']
    meta = options['meta']
    stringent_freq = options['stringent_freq']
    out_csv= [opts.out_csv[0]]
    
args = test_args()


def filter(x,args):
    stringent_variants = x.loc[ (x['ref']!=x['var']) & (x['MapQ']>args.mq) & (x['Phred']>args.phred) & (x['freq.var']>args.freq) & (x['freq.var']<args.stringent_freq) & (x['p.val']<args.pval ) & (x['Read_pos']>args.pos[0]) & (x['Read_pos']<args.pos[1])]
    not_stringent = x.loc[(x['freq.var']>args.stringent_freq)]
    stringent_ref  = x.loc[ (x['ref']==x['var']) & (x['MapQ']>args.mq) & (x['Phred']>args.phred) & (x['freq.var']>args.freq) & (x['freq.var']<args.stringent_freq) & (x['Read_pos']>args.pos[0]) & (x['Read_pos']<args.pos[1])] # same as above but no p.val
    out = stringent_variants.append(not_stringent,ignore_index=True)
    out=out.append(stringent_ref,ignore_index=True)
    out=out.drop_duplicates()
    return(out)
#print args.mq

# read in the csv and apply the subsetting function
#print args.pos
data=pd.DataFrame.from_csv(args.variants[0],index_col=False)
data.Id.apply(str)
out = filter(data,args)

if out.shape[0]>0: # if there are some variants left
    out_id=out["Id"].apply(lambda x: pd.Series(str(x).split('_')))
    out=out.assign(LAURING_ID = out_id[0])
    if len(out)==2:
        out=out.assign(dup=out_id[1])
    else:
        out=out.assign(dup=None)
    if args.run != None:
        print("Adding run label")
        out.loc[:,'run']=args.run


    if args.meta !=None:
        # split the duplicate labels if present here and then add meta data.
        print("Merging with meta data")
        #print out
        meta=pd.DataFrame.from_csv(args.meta,index_col=None)
        meta.LAURING_ID.apply(str)
        out=out.merge(meta,how="left",on="LAURING_ID")

else:
    print("There were no variants left - There may have been none to begin with")
    
print("Saving output as " + args.out_csv[0])
out.to_csv(args.out_csv[0])







