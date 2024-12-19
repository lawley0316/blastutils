# -*- coding: utf-8 -*-
from blastutils import BlastOutputFile, Reader, HSP, Hit, Record, ByHSPEvalueBitscore, MinHSPSimilarityCoverage

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


class TestSort:
    def test_hit2_empty(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref2 = Hit('ref2', 100)
        record.add(ref2)
        record.add(ref1)
        record.sort(ByHSPEvalueBitscore())
        assert record.hits[0].sseqid == 'ref1'
        assert record.hits[1].sseqid == 'ref2'

    def test_hit1_empty(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref2 = Hit('ref2', 100)
        ref2.add(HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48))
        record.add(ref1)
        record.add(ref2)
        record.sort(ByHSPEvalueBitscore())
        assert record.hits[0].sseqid == 'ref2'
        assert record.hits[1].sseqid == 'ref1'

    def test_evalue(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref1.add(HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 185.0, 8.22e-50))
        ref2 = Hit('ref2', 100)
        ref2.add(HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48))
        record.add(ref2)
        record.add(ref1)
        record.sort(ByHSPEvalueBitscore())
        assert record.hits[0].sseqid == 'ref1'
        assert record.hits[1].sseqid == 'ref2'

    def test_bitscore(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref1.add(HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 185.0, 1.06e-48))
        ref2 = Hit('ref2', 100)
        ref2.add(HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48))
        record.add(ref2)
        record.add(ref1)
        record.sort(ByHSPEvalueBitscore())
        assert record.hits[0].sseqid == 'ref1'
        assert record.hits[1].sseqid == 'ref2'


class TestBest:
    def test_empty(self):
        record = Record('seq1', 100)
        assert record.best(ByHSPEvalueBitscore()) is None

    def test_evalue(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref1.add(HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 185.0, 8.22e-50))
        ref2 = Hit('ref2', 100)
        ref2.add(HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48))
        record.add(ref2)
        record.add(ref1)
        assert record.best(ByHSPEvalueBitscore()).sseqid == 'ref1'

    def test_bitscore(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref1.add(HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 185.0, 1.06e-48))
        ref2 = Hit('ref2', 100)
        ref2.add(HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48))
        record.add(ref2)
        record.add(ref1)
        assert record.best(ByHSPEvalueBitscore()).sseqid == 'ref1'

    def test_evalue_bitscore_equal1(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref1.add(HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 182.0, 1.06e-48))
        ref2 = Hit('ref2', 100)
        ref2.add(HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48))
        record.add(ref1)
        record.add(ref2)
        assert record.best(ByHSPEvalueBitscore()).sseqid == 'ref1'

    def test_evalue_bitscore_equal2(self):
        record = Record('seq1', 100)
        ref1 = Hit('ref1', 100)
        ref1.add(HSP(1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 182.0, 1.06e-48))
        ref2 = Hit('ref2', 100)
        ref2.add(HSP(3, 100, 1625, 1722, 0, 0, 98, 100.0, 98.0, 182.0, 1.06e-48))
        record.add(ref2)
        record.add(ref1)
        assert record.best(ByHSPEvalueBitscore()).sseqid == 'ref2'


class TestFilter:
    def test(self):
        with BlastOutputFile('example.txt') as file:
            record = next(Reader(file))
        record.filter(MinHSPSimilarityCoverage(99, 99))
        assert record.count() == 1
        assert record.hits[0].sseqid == 'ref1'