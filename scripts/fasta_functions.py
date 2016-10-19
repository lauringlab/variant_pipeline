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

def Align(headers_seqs, progpath, musclegapopen=None):
    """Performs a multiple sequence alignment of two or more sequences.

    'headers_seqs' is a list of seq_record objects. The list must contain 2 objects and the reference is [first,second].

    'progpath' should specify a directory containing the alignment program executable,
        either PROBCONS or MUSCLE.  The PROBCONS executable is assumed to have
        the name "probcons" in this directory.  The MUSCLE executable is assumed to
        have the name "muscle" in this directory.


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
        raise# ValueError, 'header_seqs does not specify a list with at least two entries.'
    if not os.path.isdir(progpath):
        raise# ValueError, "Cannot find directory %s." % progpath
    exe = os.path.abspath("%s/muscle" % progpath) # the executable
    if not os.path.isfile(exe):
        raise# IOError, "Cannot find executable at %s." % exe
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
            raise#
    finally:
        os.chdir(currdir) # return to the original directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file)) # remove files from temporary directory
        os.rmdir(tempdir) # remove temporary directory
    #print headers_seqs
    if len(aligned_headers_seqs) != len(headers_seqs):
        raise# ValueError, "Did not return the correct number of aligned sequences."
    return aligned_headers_seqs # return the aligned sequences

def GetCorrespondingResidue(seqs, i):
    """Gets the corresponding residue number for two aligned sequences.

    *seqs* is a set of two aligned sequences as *(head, sequence)* 2-tuples.

    *i* is the number of a residue in sequential numbering of *seqs[0]*
    without considering any of the gaps induced by alignment, in 1, 2, ...
    numbering.

    Returns the number of the residue in sequential numbering of *seqs[1]*
    without considering any of the gaps induced by alignment in 1, 2, ...
    numbering. Returns *None* if residue *i* of *seqs[0]* aligns with a
    gap in *seqs[1]*.
    """
    assert len(seqs) == 2
    s1 = seqs[0].seq
    s2 = seqs[1].seq
    assert len(s1) == len(s2)
    assert 1 <= i <= len(s1)
    s1index = s2index = 0
    for j in range(len(s1)):
        if s1[j] != '-':
            s1index += 1
        if s2[j] != '-':
            s2index += 1
        if s1index == i:
            if s2[j] == '-':
                return None
            else:
                return s2index
def StripGapsToFirstSequence(aligned_headers_seqs):
    """Strips gaps from a reference sequence, and all corresponding alignments.
    On input, 'aligned_headers_seqs' should be a set of two or more aligned sequences,
        as would be returned by 'Align'.
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
    >>> StripGapsToFirstSequence([('s1', '-AT-A-GC'), ('s2', 'AAT-TAGC')]
    """        
    ref_seq = copy.deepcopy(aligned_headers_seqs[0].seq)
    ref_stripped=copy.deepcopy(aligned_headers_seqs[0])
    iseq_stripped = copy.deepcopy(aligned_headers_seqs[1])

    non_strip_positions = [] # positions not to strip away
    stripped_ref_seq = []
    for i in range(len(ref_seq)):
        r = ref_seq[i]
        if r != '-':
            non_strip_positions.append(i)
            stripped_ref_seq.append(r)
    
    ref_stripped.seq=Seq(''.join(stripped_ref_seq))
    ref_stripped.seq.alphabet=IUPAC.unambiguous_dna
    iseq_stripped.seq=Seq(''.join([iseq_stripped.seq[i] for i in non_strip_positions]))
    iseq_stripped.seq.alphabet=IUPAC.unambiguous_dna

    return iseq_stripped
    
def catch(func, handle=lambda e : e, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as e:
        return None
        
def mutate(sequence,variants_df):
    """ This function takes in a Seq object and a data frame with mutations with chr, pos, ref, var columns. It 
    applies the mutations and then returns a sequence containing all the mutations in the variant data frame.
    maybe I'll use lists of sequences instead of one sequence.
    """
    seq=copy.deepcopy(sequence)
    seq.seq=seq.seq.tomutable()
    
    # Get the most recent coding position 
    df=variants_df
    for index, row in df.iterrows():
        if row["ref"]!=seq.seq[int(row["pos"])-1]:
            raise ValueError("Reference base does not match the reference base in the sequence")
        seq.seq[int(row["pos"])-1]=row["var"]
    seq.seq=seq.seq.toseq()
    return seq    
    
    
def var_aa(sequence,variants_df,j):
    """ This function will take a  Seq object and a varaints data frame of mutation with chr, pos, ref, var columns.
    and the index of the row that contains the mutation in question. In the original useage, I have already fixed the
     fixed mutations. So nothing will change in the sequence when we reach that point, but that's ok. We'll still check for a difference
     This will return an updated data frame containing the varaint amino acid. It returns a dictionary with the amino Acid "AA" and the 
     class of mutation.
    """
    seq=copy.deepcopy(sequence)
    seq.seq=seq.seq.tomutable()
    try :
        seq.seq[int(variants_df.loc[j,"coding_pos"][-1])-1]=variants_df.loc[j,"var"]
    except TypeError:
        seq.seq=seq.seq
    seq.seq=seq.seq.toseq()
    protien=seq.seq.translate()
    try :
        aa=protien[int(variants_df.loc[j,"AA_pos"][-1])-1]
    except TypeError:
        aa=None
    #variants_df.loc[j,"Var_AA"].append(aa)
        
    
    if variants_df.loc[j,"Ref_AA"][-1]==aa==None:
        #variants_df.loc[j,"Class"].append("NonCoding")
        type_mut = None
    elif variants_df.loc[j,"Ref_AA"][-1]!=aa and aa=='*':
        #variants_df.loc[j,"Class"].append("Stop")
        type_mut = "Stop"
    elif variants_df.loc[j,"Ref_AA"][-1]!=aa and aa=='*':
        #variants_df.loc[j,"Class"].append("Read Through")
        type_mut = "Read Through"
    elif variants_df.loc[j,"Ref_AA"][-1]==aa:
        #variants_df.loc[j,"Class"].append("S")
        type_mut = "Syn"
    elif variants_df.loc[j,"Ref_AA"][-1]!=aa:
        #variants_df.loc[j,"Class"].append("NS")
        type_mut = "Nonsyn"
    return{"AA":aa,"Class":type_mut}
