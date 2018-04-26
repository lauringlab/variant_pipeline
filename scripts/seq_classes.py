from __future__ import division

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
    def calc_coverage(self):
        cov = []
        for loci in self.seq:
            cov.append(loci.coverage)
        return cov
     