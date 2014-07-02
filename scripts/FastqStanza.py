class FastqStanza(object):
    """ Encapsulates the four-line chunks of a fastq file. """
    @classmethod
    def parse(cls, line_str, delim="\n"):
        return FastqStanza(*line_str.rstrip().split(delim))
    
    def __init__(self, sequence_id, sequence, quality_id, quality):  
        self.sequence_id = sequence_id 
        self.sequence = sequence
        self.quality_id = quality_id
        self.quality = quality
    
    def __repr__(self):
        return "{0}\n{1}\n{2}\n{3}".format(
            self.sequence_id, self.sequence, self.quality_id, self.quality)

    def as_text(self):
        return self.__repr__() + "\n"

def stanza_generator(reader, stanza_delimiter="@"):
    #skip header
    for line in reader:
        stanza_str = line
        if stanza_str.startswith(stanza_delimiter):
            break
    count = 1       
    for line in reader:
        if count % 4 == 0:
            yield FastqStanza.parse(stanza_str)
            stanza_str = line
            count = 1
        else:
            stanza_str += line
            count += 1
    yield FastqStanza.parse(stanza_str)


def trim_to_shortest_sequence_length(fastqA, fastqB):
    a_len = len(fastqA.sequence)
    b_len = len(fastqB.sequence)
    if a_len != b_len:
        shortest_len = min(a_len, b_len)
        fastqA.sequence = fastqA.sequence[0:shortest_len]
        fastqA.quality = fastqA.quality[0:shortest_len]
        fastqB.sequence = fastqB.sequence[0:shortest_len]
        fastqB.quality = fastqB.quality[0:shortest_len]

