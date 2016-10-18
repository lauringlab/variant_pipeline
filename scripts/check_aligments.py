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


parser = argparse.ArgumentParser(description='This is a helper script that takes in a list of fasta files and a segment name and aligns them using muscle. As always I\'m borrowing from the Bloom lab\'s HA_numbering script. By defualt it prints to the screen in -clw format',usage ="")

parser.add_argument('alignerpath', metavar='alignerpath', nargs='+',
                    help='The path to the muscle executable - assuming the executable is name muscle')


parser.add_argument('-fasta', action='store',nargs="+",dest= 'fasta',
                    help='The fasta files ')
                    
parser.add_argument('-segs', action='store',nargs="+",dest= 'segs',
					help='A list of the segments you want to align. assuming the segment is named the same in all files')
					
parser.add_argument('-megaopen',action='store',dest='megaopen',default=None,type=int,
                    help='optional gap penalty for muscle. should be negative')
                    
                    
def Align(headers_seqs, progpath, musclegapopen=None):
    """Performs a multiple sequence alignment of two or more sequences.

    By default, the protein sequences are aligned using PROBCONS.  This is
        probably the most accurate alignment program.  However, it is
        slow and consumes large amounts of memory if you are aligning
        a very large number of sequences (typically if you are aligning
        more than several hundred).  In that case, you may prefer to use
        MUSCLE instead.  You can choose between the two with the 'program'
        option.  If you decide to use MUSCLE, you can also align nucleotide
        sequences with this program.

    'headers_seqs' is a list of seq_record objects. The list must contain 2 objects and the reference is [first,second].

    'progpath' should specify a directory containing the alignment program executable,
        either PROBCONS or MUSCLE.  The PROBCONS executable is assumed to have
        the name "probcons" in this directory.  The MUSCLE executable is assumed to
        have the name "muscle" in this directory.

    'program' specifies what program to use for the alignment.  By default, it is
        "PROBCONS".  If you wish to use MUSCLE instead, set it to "MUSCLE".

    'musclegapopen' sets the MUSCLE gap openining penalty to the specified
        value. By default it is None, meaning we use the MUSCLE default penalty.
        You can also set it to a number; for example -100 will lead to fewer gaps.

    This executable is used to perform a multiple sequence alignment of the proteins
        with the default settings of either PROBCONS or MUSCLE.  The returned variable is a
        new list 'aligned_headers_seqs'.  Each entry is a 2-tuple '(head, aligned_seq)'.
        'head' has the same meaning as on input (the sequence header) and
        'aligned_seq' is the aligned sequence, with gaps inserted as '-'
        as appropriate.  Therefore, all of the 'aligned_seq' entries in
        'aligned_headers_seqs' are of the same length.  The entries in 'aligned_headers_seq'
        are in the same order as in the input list 'headers_seqs'.
    """
    if not (isinstance(headers_seqs, list) and len(headers_seqs) >= 2):
        raise ValueError, 'header_seqs does not specify a list with at least two entries.'
    if not os.path.isdir(progpath):
        raise ValueError, "Cannot find directory %s." % progpath
    exe = os.path.abspath("%s/muscle" % progpath) # the executable
    if not os.path.isfile(exe):
        raise IOError, "Cannot find executable at %s." % exe
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in a temporary directory
        infile = "%s/in.fasta" % tempdir # input file
        SeqIO.write(headers_seqs, infile, "fasta") # write sequences to the input file
        if musclegapopen != None:
            p = subprocess.Popen("%s -gapopen %d -in %s -clw" % (exe, musclegapopen, infile), shell = True)#, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run MUSCLE
        else:
            p = subprocess.Popen("%s -in %s -clw" % (exe, infile), shell = True)#, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run MUSCLE
        (output, errors) = p.communicate()
    finally:
        os.chdir(currdir) # return to the original directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file)) # remove files from temporary directory
        os.rmdir(tempdir) # remove temporary directory
    
def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    'fastafile' should specify the name of a FASTA file.

    This function reads all sequences from the FASTA file.  It returns the
        list 'headers_seqs'.  This list is composed of a seq_record objects.
    """
    seqs =[]
    header = None
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq_record.seq.alphabet=IUPAC.unambiguous_dna
        seqs.append(seq_record)

    return seqs
####### Starting main script##########
args = parser.parse_args()

print args
    
for seg in args.segs:
    fa_seq=[]    
    for fastafile in args.fasta:
        seqs=ReadFASTA(fastafile)    
        for s in seqs:
            if s.id ==seg:
                s.id=s.id+str(fastafile)
                fa_seq.append(s)
    Align(fa_seq, args.alignerpath[0],args.megaopen)

    