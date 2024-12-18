# -*- coding: utf-8 -*-
from blastuitls import BlastOutputFile, Reader, HSP

EXAMPLE_FILE = 'example.txt'
ROWS = [
    ['seq1', 100, 'ref1', 100, 1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 185.0, 8.22e-50],
    ['seq1', 100, 'ref2', 1722, 3, 100, 3, 100, 0, 0, 98, 100.0, 98, 182, 1.06e-48],
    ['seq1', 100, 'ref2', 1722, 3, 100, 814, 911, 0, 0, 98, 100.0, 98, 182, 1.06e-48],
    ['seq1', 100, 'ref2', 1722, 3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48],
    ['seq2', 100, 'ref1', 100, 1, 100, 1, 100, 0, 0, 100, 100.0, 100, 185, 8.22e-50],
    ['seq2', 100, 'ref2', 1722, 3, 100, 3, 100, 0, 0, 98, 100.0, 98, 182, 1.06e-48],
    ['seq2', 100, 'ref2', 1722, 3, 100, 814, 911, 0, 0, 98, 100.0, 98, 182, 1.06e-48],
    ['seq2', 100, 'ref2', 1722, 3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48]
]


class TestBlastOutputFile:
    def test_read(self):
        with BlastOutputFile(EXAMPLE_FILE) as file:
            rows = list(file)
        assert rows == ROWS

    def test_write(self):
        with BlastOutputFile('BlastOutputFile_write_demo.txt', 'w') as file:
            for row in ROWS:
                file.write(row)


class TestReader:
    def test_total(self):
        with BlastOutputFile(EXAMPLE_FILE) as file:
            records = list(Reader(file))
        assert len(records) == 2

    def test_record0(self):
        with BlastOutputFile(EXAMPLE_FILE) as file:
            record = next(Reader(file))
        assert record.qseqid == 'seq1'
        assert record.qlen == 100
        assert record.count() == 2
        hit0, hit1 = record.hits

        # hit0
        assert hit0.sseqid == 'ref1'
        assert hit0.slen == 100
        assert hit0.count() == 1
        assert hit0.hsps[0] == HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100, 185, 8.22e-50)

        # hit1
        assert hit1.sseqid == 'ref2'
        assert hit1.slen == 1722
        assert hit1.count() == 3
        assert hit1.hsps[0] == HSP(3, 100, 3, 100, 0, 0, 98, 100.0, 98, 182, 1.06e-48)
        assert hit1.hsps[1] == HSP(3, 100, 814, 911, 0, 0, 98, 100.0, 98, 182, 1.06e-48)
        assert hit1.hsps[2] == HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98, 182, 1.06e-48)

    def test_record1(self):
        with BlastOutputFile(EXAMPLE_FILE) as file:
            record = list(Reader(file))[1]
        assert record.qseqid == 'seq2'
        assert record.qlen == 100
        assert record.count() == 2
        hit0, hit1 = record.hits

        # hit0
        assert hit0.sseqid == 'ref1'
        assert hit0.slen == 100
        assert hit0.count() == 1
        assert hit0.hsps[0] == HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100, 185, 8.22e-50)

        # hit1
        assert hit1.sseqid == 'ref2'
        assert hit1.slen == 1722
        assert hit1.count() == 3
        assert hit1.hsps[0] == HSP(3, 100, 3, 100, 0, 0, 98, 100.0, 98, 182, 1.06e-48)
        assert hit1.hsps[1] == HSP(3, 100, 814, 911, 0, 0, 98, 100.0, 98, 182, 1.06e-48)
        assert hit1.hsps[2] == HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98, 182, 1.06e-48)