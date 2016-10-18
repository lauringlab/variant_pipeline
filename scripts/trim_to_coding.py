#!/Users/jt/.virtualenvs/sci-py2.7/bin/python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import argparse
import glob
import os
import tempfile
import subprocess


""" This script takes in 2 fasta files and trimms the sequences in 1 to that of the reference.
    I am using MUSCLE to align the sequences, and am drawing heavily (in many cases verbatim) from the HA_number script
    deleveped by Jesse Bloom at https://github.com/jbloomlab/HA_numbering/blob/master/HA_numbering.py
"""


parser = argparse.ArgumentParser(description='This script takes in a fasta test fasta file and trims it to match the regions found a reference fasta file.\n I am using it to trim whole genomes to just the coding regions, but I suppose it could have other uses. Currently it relies on MUSCLE. The segment names in the sample file must match those in the reference.')
parser.add_argument('aligner_path', metavar='aligner_path', nargs='+',
                    help='The path to the muscle executable - assuming the executable is name muscle')

parser.add_argument('in_fa', metavar='in_fa', nargs='+',
                    help='The input (sample) fa')

parser.add_argument('ref_fa', metavar='ref', nargs='+',
                    help='The reference fasta to which the sequences will be trimmed.')

parser.add_argument('-out_fa',action='store',dest='out_fa',default=None,
                    help='optional output the trimmed fasta file')

parser.add_argument('-csv',action='store',dest='csv',default=None,
                    help='optional output - a csv file recording the number of bp trimmed off the 5\' and 3\' ends')

args = parser.parse_args()

csv=args.csv

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
        outfile = "%s/out.fasta" % tempdir # output file
        SeqIO.write(headers_seqs, infile, "fasta") # write sequences to the input file
        if musclegapopen != None:
            p = subprocess.Popen("%s -gapopen %d -in %s -out %s" % (exe, musclegapopen, infile, outfile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run MUSCLE
        else:
            p = subprocess.Popen("%s -in %s -out %s" % (exe, infile, outfile), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run MUSCLE
        (output, errors) = p.communicate()
        try:
            aligned_headers_seqs = ReadFASTA(outfile)
        except:
            sys.stderr.write("Error getting alignment output, error of %s" % errors)
            raise
    finally:
        os.chdir(currdir) # return to the original directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file)) # remove files from temporary directory
        os.rmdir(tempdir) # remove temporary directory
    if len(aligned_headers_seqs) != len(headers_seqs):
        raise ValueError, "Did not return the correct number of aligned sequences."
    # put the aligned sequences in the same order as the input sequences
    # n = len(aligned_headers_seqs[0][1]) # length of aligned sequences
    # d = dict(aligned_headers_seqs)
    # aligned_headers_seqs = []
    # for (head, seq) in headers_seqs:
    #     try:
    #         alignedseq = d[head]
    #     except KeyError:
    #         raise ValueError("After alignment, the following header is missing: %s" % head)
    #     if len(alignedseq) != n:
    #         open('errors.temp', 'w').write(errors)
    #         raise ValueError("Aligned sequence %s is not of length %d: if you are using MUSCLE, you may be running out of memory.  Errors have been written to errors.temp." % (alignedseq, n))
    #     if len(seq) > n:
    #         open('errors.temp', 'w').write(errors)
    #         raise ValueError("Unaligned seq %s is longer than aligned length of %d: if you are using MUSCLE, you many be running out of memory.  Errors have been written to errors.temp." % (seq, n))
    #     aligned_headers_seqs.append((head, alignedseq))
    #print(aligned_headers_seqs)
    return aligned_headers_seqs # return the aligned sequences


def trim(aligned_headers_seqs):
    """trimms excess 5' and 3' regions from sample seq.


    The first sequence in this alignment is taken to correspond to the reference sequence.
        The returned variable is a list similar to aligned_headers_seqs, but with
        all positions corresponding to gaps in this reference sequence stripped away.
        All gaps ('-') characters are removed from this reference sequence.  In addition,
        in all other aligned sequences in 'aligned_headers_seqs', every character at
        the same position as a gap in the reference sequence is removed.  Therefore,
        at the end of this procedure, all of the alignments have the same length
        as the reference sequence with its gaps stripped away.  The headers are
        unchanged.  The order of sequences in this stripped alignment is also
        unchanged.

        I have test the results and the output is correct
    """
    if not (isinstance(aligned_headers_seqs, list) and len(aligned_headers_seqs) >= 2):
        raise ValueError, "Input does not specify at least two aligned sequences."
    ref_seq = aligned_headers_seqs[0].seq# str yields the sequence
    #print(ref_seq)
    # Getting the positions to strip from the start
    go=True
    i=0
    start_excess=0
    while (go==True):
        if (ref_seq[i]=='-'):
            start_excess=i # strip 0 to i
        else:
            go=False
        i=i+1
    # Getting the posisiton to remove from the end
    start_excess=start_excess+1 # slicing is inclusive on this end
    end=True
    i=len(ref_seq)-1
    end_excess=i
    print(i)
    while (end==True):
        if (ref_seq[i]=='-'):
            end_excess=i # strip 0 to i
        else:
            end=False
        i=i-1

    print "%s bases taken off the 5' end" % str(start_excess)
    print "%s bases taken off the 3' end " % str(len(ref_seq)-1-end_excess)



    samp_seq=aligned_headers_seqs[1]
    samp_seq.seq=samp_seq.seq[start_excess:end_excess]

    return([samp_seq,start_excess,end_excess+1]) # In a 1 base system (like R) The start will be the last base to not be exclued on the 5' and end is the last base off the end to be included.



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


def main(): # The positions will be given as base 0 and adjusted to match the convention (base 1) in the funciton
    """Main body of script."""
    print "\nBeginning execution trimming script."

    # parse arguments
    args = parser.parse_args()
    alignerpath = args.aligner_path[0]
    if not os.path.isdir(alignerpath):
        raise IOError,"The directory of %s specified by musclepath does not exist." % (alignerpath)
    prog = 'MUSCLE'

    sample=ReadFASTA(args.in_fa[0])
    ref=ReadFASTA(args.ref_fa[0])

    samp_seqname=[]
    ref_seqname=[]
    for seq in sample:
        samp_seqname.append(seq.id)
    for seq in ref:
        ref_seqname.append(seq.id)
    # make alignments
    print("Making %s alignments..." % prog)
    align_ref = []
    align_samp=[]
    for seqname in samp_seqname:
        #print("Aligning %s" % seqname)
        sample_seq=sample[samp_seqname.index(seqname)]
        try:
            ref_seq=ref[ref_seqname.index(seqname)]
        except ValueError:
            raise ValueError, " Segement %s was not found in the reference sequence" % seqname

        alignments=Align([ref_seq, sample_seq], alignerpath)
        align_ref.append(alignments[0])
        align_samp.append(alignments[1])
    print("Trimming...\n")
    trimmed=[]
    segs=[]
    off_5=[]
    off_3=[]
    for i in range(0,len(align_samp)):
        print "Trimming %s" % align_samp[i].id
        trimmed_out=trim([align_ref[i],align_samp[i]])
        trimmed.append(trimmed_out[0])
        segs.append(align_samp[i].id)
        off_5.append(trimmed_out[1])
        off_3.append(trimmed_out[2])





    if(csv==None):
        print "writing output to %s"  % args.out_fa
        SeqIO.write(trimmed, args.out_fa, "fasta")
    else:
        print "writing csv file to %s" % csv
        with open(csv,'w') as out_file:
           out_file.write("chr,off.5,off.3\n")
           for i in range(0,len(off_5)) :
               out_file.write(str(segs[i])+","+str(off_5[i])+","+str(off_3[i])+'\n')


main()
