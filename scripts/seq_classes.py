from __future__ import division
import pysam
class locus(object):
    """ A base which has the following characteristics 
        Atrributes :
            Chr : chromosome
            pos : position on chromosome
            counts : a dictionary of { A : the number of As
            T : the number of Ts
            C : the number of Cs
            G : the number of Gs
            - : the number of of deletions    
            }
            add consensus.
        methods : 
            update : update counts
            calc_freqs :  calculate the frequency of each base
            consensus : caculate the consensus sequence as this position 
    """
    def __init__(self,chr,pos):
        """
        return a base object with chr and pos and starting counts of 0
        Posisition is base 1 like most biologists are use to.
        """
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
    def consensus(self,cutoff=None):
        self.calc_freqs()
        v=list(self.freqs.values())
        k=list(self.freqs.keys())
        # check for cutoff method
        if cutoff ==None:
            # Return the most common base
            return k[v.index(max(v))]
        else:
            # return the base that is above the cutoff
            if type(cutoff)!=float or (cutoff>1 or cutoff<0.5):
                raise ValueError('cutoff must be a float in [0.5,1.0] or nothing')
            indexes = [index for index, value in enumerate(v) if value > cutoff]
            if len(indexes)==1:
                consensus = k[indexes[0]]
            else:
                consensus = "N"
            return consensus

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
        if len(self.seq)>0 and loci.pos!=max([x.pos for x in self.seq])+1:
            raise ValueError('The position of the loci does not match current segment length')
        self.seq.append(loci)
    def consensus(self,cutoff=None):
        seg_consensus = ""
        for loci in self.seq:
            seg_consensus=seg_consensus+loci.consensus(cutoff)
        return seg_consensus
    def locus(self,pos):
        loci = [x for x in self.seq if x.pos==pos]
        return loci[0]
    def calc_coverage(self):
        cov = []
        for loci in self.seq:
            cov.append(loci.coverage)
        return cov

# reference object

def tally(bamfile,chr,start,stop,maxDepth = 1000,phred = 30):
    """
        This function takes in a pysam alignment file and tallies base calls 
        across a defined region. Bases must pass a phred score cutoff.
    """
    seg = segment(chr)
    pileup = bamfile.pileup(chr,start,stop,stepper='all',truncate=True,max_depth = maxDepth)
#        pileup = bamfile.pileup(chr,pos,pos+1,stepper='all',truncate=True,max_depth = maxDepth)
    for pileupColumn in pileup:
        l = locus(chr,pos=pileupColumn.pos+1) # loci positions are base 1
        for pileupRead in pileupColumn.pileups:
            if not pileupRead.is_del:
                if pileupRead.is_refskip:
                    l.update("-")
                else:
                    # query_alignment_* removes solf clipped bases
                    called_base= pileupRead.alignment.query_alignment_sequence[pileupRead.query_position] 
                    called_phred= pileupRead.alignment.query_alignment_qualities[pileupRead.query_position]
                    if called_phred>=phred: 
                        l.update(called_base)
        
        seg.append_loci(l)
    return seg
    