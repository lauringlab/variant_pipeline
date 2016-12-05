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

# parser.add_argument('-mapping', action='store',dest= 'mq',type=float,default=0,
#                     help='The average mapping quality cut off')
# 
# parser.add_argument('-phred', action='store',dest= 'phred',type=float,default=0,
#                     help='The average phred quality cut off')
# 
# parser.add_argument('-freq', action='store',dest= 'freq',type=float,default=0,
#                     help='The frequency cut off')
# 
# parser.add_argument('-pos', action='store',type=int,nargs="+",dest= 'pos',default=[0,250],
#                     help='[min max] mean positions on the read')
# 
# parser.add_argument('-infer', action='store_true',default='False',dest= 'infer',
#                     help='Boolean switch to infer minor variants at positions where the minor variant may be the plasmid\'s consensus')

# parser.add_argument('-p', action='store',type=float,default=0.01,dest= 'pval',
#                     help='the p value threshold : default = 0.01')
# 
# parser.add_argument('-run',action='store',type=str,dest='run',
#                     help = 'The name of the run you would like to add to the data. If you\'d like')
# 
# parser.add_argument('-meta',action='store',dest = 'meta',
# 			help = 'The meta data file')


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
    infer = options['infer']
    out_csv= [opts.out_csv[0]]
    
args = test_args()

def infer(x):
    x=x[:1]
    x["var"].update(x["ref"])
    x["freq.var"].update(1-x["total_freq"])
    x["mutation"].update(x["chr"] +"_"+ x["ref"] +x["pos"].map(str)+x["var"])
    return(x)    


def infer_all(data,low,high):
    data=data.assign(id_pos = data['chr']+ data['pos'].map(str)+"-"+data["Id"].map(str))
    data=data.assign(exp_freq = (data['n.tst.bw']+data['n.tst.fw'])/(data['cov.tst.bw']+data['cov.tst.fw'])) # This is the expected frequency if there was nothing in the control
    total_freq=data.groupby("id_pos")['exp_freq'].apply(np.sum)
    total_freq=total_freq.to_frame("total_freq")
    total_freq.loc[total_freq.index,"id_pos"]=total_freq.index
    #print total_freq
    data=data.merge(total_freq,how='left',on='id_pos')
    #print(data.columns)
  
    to_infer=data[(data["total_freq"]>low) & (data["total_freq"]<high)]
  
    infered=to_infer.groupby("id_pos").apply(infer)
    #print(infered)
    data=data.append(infered,ignore_index=True)
    #print(data)
    data=data.drop('id_pos', axis=1)
    return(data)


def filter(x,args):
    out = x.loc[(x['MapQ']>args.mq) & (x['Phred']>args.phred) & (x['freq.var']>args.freq) & (x['p.val']<args.pval) & (x['Read_pos']>args.pos[0]) & (x['Read_pos']<args.pos[1])]
    return(out)
#print args.mq

# read in the csv and apply the subsetting function
print args.pos
data=pd.DataFrame.from_csv(args.variants[0],index_col=False)
data.Id.apply(str)
out = filter(data,args)

if out.shape[0]>0:
    out_id=out["Id"].apply(lambda x: pd.Series(str(x).split('_')))
    out=out.assign(LAURING_ID = out_id[0])
    if len(out)==2:
        out=out.assign(dup=out_id[1])
    else:
        out=out.assign(dup=None)

    if args.infer == True:
        print "Infering reciprocal variants"
        out=infer_all(out,0.01,0.99)

    if args.run != None:
        print "Adding run label"
        out.loc[:,'run']=args.run


    if args.meta !=None:
        # split the duplicate labels if present here and then add meta data.
        print "Merging with meta data"
        #print out
        meta=pd.DataFrame.from_csv(args.meta,index_col=False)
        meta.LAURING_ID.apply(str)
        out=out.merge(meta,how="left",on="LAURING_ID")

else:
    print "There were no variants left - There amy have been none to begin with"
    
print "Saving output as " + args.out_csv[0]
out.to_csv(args.out_csv[0])







