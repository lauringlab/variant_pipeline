import itertools
import os
import sys

import FastqStanza

def trim(source_left, source_right, dest_left, dest_right):
    fastq_pairs = itertools.izip(
        FastqStanza.stanza_generator(source_left), FastqStanza.stanza_generator(source_right))
    for (left_stanza, right_stanza) in fastq_pairs:
        FastqStanza.trim_to_shortest_sequence_length(left_stanza, right_stanza)
        dest_left.write(left_stanza.as_text())
        dest_right.write(right_stanza.as_text())
