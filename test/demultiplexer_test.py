import unittest

from scripts.demultiplexer import demultiplex, MismatchedSequenceIdError, UndefinedBarcodeError

class DemultiplexerTest(unittest.TestCase):
    def test_demultiplex_singleton_identity(self):
        left_stanza = "@CAT:DUPE_3:ID_1047:FLAG_1 1\nACGTACGT\n+\nAAAAAAAA"
        right_stanza = "@CAT:DUPE_3:ID_1047:FLAG_1 2\nACGTACGT\n+\nAAAAAAAA"
        source_left = MockReader(left_stanza)
        source_right = MockReader(right_stanza)
        dest_left = MockWriter()
        dest_right = MockWriter()
        barcode_files = {'CAT' : [dest_left, dest_right] }

        demultiplex(source_left, source_right, barcode_files)

        self.assertEqual(left_stanza.splitlines(), dest_left.lines())
        self.assertEqual(right_stanza.splitlines(), dest_right.lines())


    def test_demultiplex(self):
        left_stanza1 = "@AAACCC:DUPE_3:ID_1047:FLAG_1 1\nACGTACGT\n+\nAAAAAAAA"
        right_stanza1 = "@AAACCC:DUPE_3:ID_1047:FLAG_1 2\nACGTACGT\n+\nAAAAAAAA"
        left_stanza2 = "@GGGTTT:DUPE_3:ID_1047:FLAG_1 1\nACGTACGT\n+\nAAAAAAAA"
        right_stanza2 = "@GGGTTT:DUPE_3:ID_1047:FLAG_1 2\nACGTACGT\n+\nAAAAAAAA"
        source_left = MockReader(left_stanza1 + "\n" + left_stanza2)
        source_right = MockReader(right_stanza1 + "\n" + right_stanza2)
        left_barcodeA = MockWriter()
        right_barcodeA = MockWriter()
        left_barcodeB = MockWriter()
        right_barcodeB = MockWriter()
        barcode_files = { \
            'AAACCC' : [left_barcodeA, right_barcodeA], \
            'GGGTTT' : [left_barcodeB, right_barcodeB] }

        demultiplex(source_left, source_right, barcode_files)

        self.assertEqual(left_stanza1.splitlines(), left_barcodeA.lines())
        self.assertEqual(right_stanza1.splitlines(), right_barcodeA.lines())
        self.assertEqual(left_stanza2.splitlines(), left_barcodeB.lines())
        self.assertEqual(right_stanza2.splitlines(), right_barcodeB.lines())


    def test_demultiplex_fails_on_unrecognized_barcode(self):
        left_stanza = "@CAT:DUPE_3:ID_1047:FLAG_1 1\nACGTACGT\n+\nAAAAAAAA"
        right_stanza = "@CAT:DUPE_3:ID_1047:FLAG_1 2\nACGTACGT\n+\nAAAAAAAA"
        source_left = MockReader(left_stanza)
        source_right = MockReader(right_stanza)
        barcode_files = {'TTT' : [MockWriter(), MockWriter()] }

        self.assertRaises(UndefinedBarcodeError, demultiplex, source_left, source_right, barcode_files)


    def test_demultiplex_fails_on_mismatched_headers(self):
        left_stanza = "@CAT:DUPE_3:ID_1047:FLAG_1 1\nACGTACGT\n+\nAAAAAAAA"
        right_stanza = "@GTA:DUPE_3:ID_1047:FLAG_1 2\nACGTACGT\n+\nAAAAAAAA"
        source_left = MockReader(left_stanza)
        source_right = MockReader(right_stanza)
        barcode_files = {'CAT' : [MockWriter(), MockWriter()] }

        self.assertRaises(MismatchedSequenceIdError, demultiplex, source_left, source_right, barcode_files)


class MockReader():
    def __init__(self, content):
        lines = [line + "\n" for line in content.split("\n") if line != ""]
        self._iter = lines.__iter__()
        self.wasClosed = False

    def __iter__(self):
        return self._iter

    def next(self):
        self._iter.next()

    def close(self):
        self.wasClosed=True


class MockWriter():
    def __init__(self):
        self._content = []
        self.wasClosed = False

    def write(self, content):
        self._content.append(content)
        
    def lines(self):
        return "".join(self._content).splitlines()

    def close(self):
        self.wasClosed = True
