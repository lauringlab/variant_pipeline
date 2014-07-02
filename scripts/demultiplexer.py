import itertools
import os
import sys

import FastqStanza

class DemultiplexError(Exception):
    """Base class for exceptions in this module."""
    pass

class UndefinedBarcodeError(DemultiplexError):
    def __init__(self, stanza):
        super(UndefinedBarcodeError, self).__init__()
        self.sequence_id = stanza.sequence_id
    
    def __str__(self):
        return repr("sequence_id [{0}] is not in recognized barcodes." \
            .format(self.sequence_id))


class MismatchedSequenceIdError(DemultiplexError):
    def __init__(self, left_stanza, right_stanza):
        super(MismatchedSequenceIdError, self).__init__()
        self.left_sequence_id = left_stanza.sequence_id
        self.right_sequence_id = right_stanza.sequence_id

    def __str__(self):
        return repr("left/right sequence_ids do not match [{0}] != [{1}].". \
            format(self.left_sequence_id), self.right_sequence_id)


def demultiplex(left_fastq, right_fastq, barcode_files):
    fastq_pairs = itertools.izip(
        FastqStanza.stanza_generator(left_fastq), FastqStanza.stanza_generator(right_fastq))
    for (left_stanza, right_stanza) in fastq_pairs:
        assert_sequence_ids_match(left_stanza, right_stanza)
        barcode = left_stanza.sequence_id.split(':')[0][1:]
        try:
            (left_out, right_out) = barcode_files[barcode]
        except KeyError:
            raise UndefinedBarcodeError(left_stanza)
        left_out.write(left_stanza.as_text())
        right_out.write(right_stanza.as_text())

def assert_sequence_ids_match(left_stanza, right_stanza):
    if left_stanza.sequence_id.split(' ')[0]!=right_stanza.sequence_id.split(' ')[0]:
        raise MismatchedSequenceIdError(left_stanza, right_stanza)

