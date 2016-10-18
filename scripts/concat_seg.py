#!/Users/jt/.virtualenvs/sci-py2.7/bin/python

from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import glob
import pandas as pd
import os

parser = argparse.ArgumentParser(description='This script is designed to concatenate the designated segment from a director of sample fasta files')


parser.add_argument('dir', metavar='dir', nargs='+',
                    help='The directory that contains the fasta files')

parser.add_argument('seg', metavar='seg', nargs='+',
                    help='The segment to be taken from each file' )
parser.add_argument('key', metavar='key', nargs='+',
                    help='The Key to be used to add metadata to each sequence' )

parser.add_argument('out_fa', metavar='out_fa', nargs='+',
                    help='The output file')


args = parser.parse_args()




def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    'fastafile' should specify the name of a FASTA file.

    This function reads all sequences from the FASTA file.  It returns the
        list 'headers_seqs'.  This list is composed of a seq_record objects.
    """
    seqs =[]
    header = None
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seqs.append(seq_record)

    return seqs



meta_df=pd.read_csv(args.key[0])
lauring_id=list(meta_df["Id"])

print lauring_id


selected_seg=[]
for fa in glob.glob(args.dir[0]+"/*.fa"):
    seqs=ReadFASTA(fa)
    for seq in seqs:
        if (seq.id==args.seg[0]):
            Id=os.path.basename(fa).split('.')[0]
            try:
                if int(Id) in lauring_id:
                    i=lauring_id.index(int(Id))
                    intervention=list(meta_df["Intervention"])[i]
                    if args.seg[0]=="HA":
                        geom=int(list(meta_df["HAI.geo"])[i])
                    if args.seg[0]=="NR":
                        geom=int(list(meta_df["NAI.geo"])[i])
                    date=list(meta_df["collection_date"])[i]
                    print Id
                else:
                    print "didn't find" + Id
                    intervention="NA"
                    geom="NA"
                    date="NA"
            except ValueError:
                    print "value error for " + Id
                    intervention="NA"
                    geom="NA"
                    date="NA"

            seq.id=str(Id)+"_"+str(intervention)+"_"+str(geom)+"_"+str(date)
            seq.description=""
            selected_seg.append(seq)
SeqIO.write(selected_seg, args.out_fa[0], "fasta")
