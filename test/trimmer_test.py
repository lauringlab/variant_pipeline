import unittest

from scripts.trimmer import trim

class TrimmerTest(unittest.TestCase):
    def test_trim_passthroughIfLengthsMatch(self):
        left_stanza = "@id.1\nACGTACGT\n+\nQQQQQQQQ"
        right_stanza = "@id.2\nTGCATGCA\n+\nQQQQQQQQ"
        source_left = MockReader(left_stanza)
        source_right = MockReader(right_stanza)
        dest_left = MockWriter()
        dest_right = MockWriter()

        trim(source_left, source_right, dest_left, dest_right)

        self.assertEqual(left_stanza.splitlines(), dest_left.lines())
        self.assertEqual(right_stanza.splitlines(), dest_right.lines())


    def test_trim_trimsToMinLength(self):
        left_stanza = "@id.1\nACGT\n+\nQQQQ"
        right_stanza = "@id.2\nTGCATGCA\n+\nQQQQQQQQ"
        source_left = MockReader(left_stanza)
        source_right = MockReader(right_stanza)
        dest_left = MockWriter()
        dest_right = MockWriter()

        trim(source_left, source_right, dest_left, dest_right)

        self.assertEqual(left_stanza.splitlines(), dest_left.lines())
        self.assertEqual(['@id.2', 'TGCA', '+', 'QQQQ'], dest_right.lines())


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
