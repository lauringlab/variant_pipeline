from __future__ import division
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from trim_to_coding import trim_to_regions


def checkORF(seq):
    """
    This function varifies the input is an ORF by checking it has
    i) a start codon
    ii) is made of complete codons mod(lenth(seg),3)=0 taking gaps into account
    iii) Has a stop codon at the end and not before.
    """
    protien = seq.translate()
   # check for ATG at start
    if seq.find('ATG') !=0:
        if seq.find('ATG') ==-1:
            raise ValueError("No start codon")
        if seq.find('ATG')>0:
            raise ValueError("start codon found at position "+ str(seq.find('ATG')))
    
    elif len(seq)% 3!=0:
        raise ValueError("The sequence is not multiple of 3") 
       
    elif protien.find('*')  != len(protien)-1:
        if protien.find('*') ==-1:
            raise ValueError("No stop codon in ORF")
        else:
            raise ValueError("Internal stop codon found at position "+ str(protien.find('*')))
    
    else:
        return(True)

def classify(sequence,codingRegion,pos,nucleotide):
       
        """
        seqeunce is a seqRecord
        codingRegion is a diction in the form 
        {
            "name": "NS1",
            "regions": [
                {
                    "start": 26,
                    "stop": 719
                }
            ]
        }
        The ouput is a dictionary added to the mutationalClass list
        it is of the form 
        {
            ORF: The name of the ORF,
            codonPos: the position in the codon either [0,1,2],
            codingPos: the nucleotide position in the ORF
            aminoAcidPos: The position of the amino acid in the polypetide,
            consensusAA: The consensus amino acid,
            varAA: The variant amino acid,
            classification: Nonsynonymous,Synonymous,indel, stop,
        }
        """
        # Get the coding sequence
        # Get the new position in the coding sequence
        i = 0
        # Has to catch the case where it's outside these regions
        outsideORF=True
        for seg in codingRegion["regions"]:
            if pos >= seg["start"] and pos<seg["stop"]:
                outsideORF=False
                break
            else:
                i+=1
        if outsideORF:
            return({
                "ORF": codingRegion["name"],
                "codonPos": None,
                "codingPos": None,
                "aminoAcidPos": None,
                "consensusAA": None,
                "varAA": None,
                "classification": "Noncoding"

            })
        
        codingSequence=""
        for seg in codingRegion["regions"]:
            codingSequence=codingSequence+sequence.seq[seg["start"]:seg["stop"]]
        consensusSequence = Seq(codingSequence,generic_dna)
        checkORF(consensusSequence)
        
        posInSeg = pos-codingRegion["regions"][i]["start"]
        otherRegions = codingRegion["regions"][:i]
        adjustment = 0
        if len(otherRegions)>0:
            for seg in otherRegions:
                adjustment+=seg["stop"]-seg["start"]
        
        codingPos = posInSeg+adjustment
        codonPos = codingPos % 3
        aminoAcidPos = codingPos // 3 

        consensusProtien = consensusSequence.translate()
        consensusAA = consensusProtien[aminoAcidPos]

        if nucleotide=="-":
            return({
                "ORF": codingRegion["name"],
                "codonPos": codonPos,
                "codingPos": codingPos,
                "aminoAcidPos": aminoAcidPos,
                "consensusAA": consensusAA,
                "varAA": None,
                "classification": "Indel"
            })

        mutantCodingSequence = codingSequence        
        mutantCodingSequence = mutantCodingSequence[:codingPos]+ str(nucleotide) + mutantCodingSequence[codingPos+1:]
        mutantSequence = Seq(mutantCodingSequence,generic_dna)
        mutantProtein = mutantSequence.translate()

        varAA = mutantProtein[aminoAcidPos]

        if varAA==consensusAA:
            classification = "Synonymous"
        elif varAA=="*":
            classification = "Stop"
        elif consensusAA=="*" and varAA!=consensusAA:
            classification= "Readthrough"
        elif varAA!=consensusAA:
            classification="Nonsynonymous"

        return({
                "ORF": codingRegion["name"],
                "codonPos": codonPos,
                "codingPos": codingPos,
                "aminoAcidPos": aminoAcidPos,
                "consensusAA": consensusAA,
                "varAA": varAA,
                "classification": classification
                 })



class allele(object):
    """
    The allele present at a loci and their accompanying data
    """
    def __init__(self,nucleotide):
        self.nucleotide=nucleotide
        self.count = 0
        self.freq = 0
        self.mutationalClass =[] 
    def classifyVar(self,sequence,codingRegion,pos):
        self.mutationalClass.append(classify(sequence,codingRegion,pos,self.nucleotide))
    def reprJSON(self): # https://stackoverflow.com/questions/5160077/encoding-nested-python-object-in-json
        d = dict()
        for a, v in self.__dict__.items():
            if (hasattr(v, "reprJSON")):
                d[a] = v.reprJSON()
            else:
                d[a] = v
        return d

    


    


class locus(object):
    """ A base which has the following characteristics 
        Atrributes :
            Chr : chromosome
            pos : position on chromosome base 0
            counts : a dictionary of { A : the number of As
            T : the number of Ts
            C : the number of Cs
            G : the number of Gs
            - : the number of of deletions    
            }
            add consensus.
            add concat_pos
        methods : 
            update : update counts
            calc_freqs :  calculate the frequency of each base
            consensus : caculate the consensus sequence as this position 
    """
    def __init__(self,pos):
        """
        return a base object with chr and pos and starting counts of 0
        Posisition is base 0 because this is python and I want to keep everything 
        easy to remember. If it came from python it is base 0.
        """
        self.pos = pos
        self.alleles = {'A':allele("A"),
                        'T':allele("T"),
                        'C':allele("C"),
                        'G':allele("G"),
                        '-':allele("-")}
        self.coverage = 0
        self.concat_pos = None
        self.consensus = ''

    def update(self,base):
        self.coverage = self.coverage+1

        self.alleles[base].count += 1
        self.calc_freqs()
        self.consensus = self.calc_consensus()
    def calc_freqs(self):
        for base in self.alleles.keys():
            self.alleles[base].freq = self.alleles[base].count/self.coverage
            
    def calc_consensus(self,cutoff=None):
        self.calc_freqs()
        v=[y.freq for y  in self.alleles.values()] #values
        k=list(self.alleles.keys()) # key
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

    def reprJSON(self): # https://stackoverflow.com/questions/5160077/encoding-nested-python-object-in-json
        d = dict()
        for a, v in self.__dict__.items():
            if (hasattr(v, "reprJSON")):
                d[a] = v.reprJSON()
            elif a=="alleles":
                d["alleles"]={}
                print v
                for nt, alle in v.items():
                    d["alleles"].update({nt: alle.reprJSON()})
            else:
                d[a] = v
        return d

class segment(object):
    """ A sequence like object made up of locus objects
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
        # if loci.chr!=self.chr:
        #     raise ValueError('The loci chr does not match the segement chr')
        if len(self.seq)>0 and loci.pos!=max([x.pos for x in self.seq])+1:
            raise ValueError('The position of the loci does not match current segment length')
        self.seq.append(loci)
    def consensus(self,cutoff=None):
        seg_consensus = ""
        for loci in self.seq:
            seg_consensus=seg_consensus+loci.calc_consensus(cutoff)
        return seg_consensus
    def locus(self,pos):
        loci = [x for x in self.seq if x.pos==pos]
        return loci[0]
    def calc_coverage(self):
        cov = []
        for loci in self.seq:
            cov.append(loci.coverage)
        return cov
    
    def reprJSON(self): # https://stackoverflow.com/questions/5160077/encoding-nested-python-object-in-json
        d = dict()
        for a, v in self.__dict__.items():
            if hasattr(v, "reprJSON"):
                d[a] = v.reprJSON()
            elif a=="seq":
                d["seq"]=[]
                for l in v:
                    d["seq"].append(l.reprJSON())
            else:
                d[a] = v
        return d

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
        l = locus(pos=pileupColumn.pos) # loci positions are base 0
        for pileupRead in pileupColumn.pileups:
            if not pileupRead.is_del:
                if pileupRead.is_refskip:
                    l.update("-")
                else:
				# query_alignment_* removes solf clipped bases, this confirms we don't look at positions beyond this clipping

                    if pileupRead.query_position >= pileupRead.alignment.query_alignment_start and pileupRead.query_position < pileupRead.alignment.query_alignment_end:
                    	try:
                    	    called_base= pileupRead.alignment.query_alignment_sequence[pileupRead.query_position-pileupRead.alignment.query_alignment_start] 
                    	    called_phred= pileupRead.alignment.query_alignment_qualities[pileupRead.query_position-pileupRead.alignment.query_alignment_start]
                        except IndexError:
                            print " chr : %s [%d - %d] query_pos : %d pos: %d" %(chr,pileupRead.alignment.query_alignment_start,pileupRead.alignment.query_alignment_end,pileupRead.query_position,pileupColumn.pos)
                            print pileupRead.alignment.query_sequence
                            print "-"*pileupRead.alignment.query_alignment_start + pileupRead.alignment.query_alignment_sequence
                            print pileupRead.alignment.is_reverse
                        if called_phred>=phred: 
                            l.update(called_base)
        
        seg.append_loci(l)
    return seg
    