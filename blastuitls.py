# -*- coding: utf-8 -*-
from pathlib import Path


class HSP:
    """High-scoring Segment Pair.

    Args:
        qstart: <int>, Start of alignment in query
        qend: <int>, End of alignment in query
        sstart: <int>, Start of alignment in subject
        send: <int>, End of alignment in subject
        mismatch: <int>, Number of mismatches
        gapopen: <int>, Number of gap openings
        length: <int>, Alignment length
        pident: <float>, Percentage of identical matches
        qcovhsp: <float>, Query Coverage Per HSP
        bitscore: <float>, Bit score
        evalue: <float>, Expect value
    """
    def __init__(self, qstart, qend, sstart, send, mismatch, gapopen, length, pident, qcovhsp, bitscore, evalue):
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.length = length
        self.pident = pident
        self.qcovhsp = qcovhsp
        self.bitscore = bitscore
        self.evalue = evalue

    def __str__(self):
        return (
            f'{self.qstart}\t{self.qend}\t{self.sstart}\t{self.send}\t{self.mismatch}\t'
            f'{self.gapopen}\t{self.length}\t{self.pident}\t{self.qcovhsp}\t{self.bitscore}\t{self.evalue}'
        )

    def __repr__(self):
        return f'<HSP: ({self.qstart}, {self.qend}) on ({self.sstart}, {self.send})>'

    def __eq__(self, other):
        if not isinstance(other, HSP):
            return False
        return (
            self.qstart == other.qstart) and (self.qend == other.qend) and (
            self.sstart == other.sstart) and (self.send == other.send) and (
            self.mismatch == other.mismatch) and (self.gapopen == other.gapopen) and (
            self.length == other.length) and (self.pident == other.pident) and (
            self.qcovhsp == other.qcovhsp) and (self.bitscore == other.bitscore) and (
            self.evalue == other.evalue
        )


class Hit:
    """Each Hit object represents one BLAST hit of the query against a subject.

    Args:
        sseqid: <str>, Subject Seq-id
        slen: <int>, Subject sequence length
    """
    def __init__(self, sseqid, slen):
        self.sseqid = sseqid
        self.slen = slen
        self.hsps = []

    def add(self, hsp):
        self.hsps.append(hsp)

    def count(self):
        return len(self.hsps)

    def is_empty(self):
        return self.count() == 0

    def create(self, qstart, qend, sstart, send, mismatch, gapopen, length, pident, qcovhsp, bitscore, evalue):
        hsp = HSP(qstart, qend, sstart, send, mismatch, gapopen, length, pident, qcovhsp, bitscore, evalue)
        self.add(hsp)
        return hsp

    def __str__(self):
        if self.is_empty():
            return ''
        return '\n'.join(f'{self.sseqid}\t{self.slen}\t{hsp}' for hsp in self.hsps)

    def __repr__(self):
        return f'<Hit for {self.sseqid}>'


class Record:
    """A Record object stores the information provided by BLAST for a single query.

    Args:
        qseqid: <str>, Query Seq-id
        qlen: <int>, Query sequence length
    """
    def __init__(self, qseqid, qlen):
        self.qseqid = qseqid
        self.qlen = qlen
        self.hits = []

    def add(self, hit):
        self.hits.append(hit)

    def count(self):
        return len(self.hits)

    def is_empty(self):
        return self.count() == 0

    def _is_last(self, sseqid):
        if self.is_empty():
            return False
        return self.hits[-1].sseqid == sseqid

    def create(self, sseqid, slen, qstart, qend, sstart, send, mismatch, gapopen, length, pident, qcovhsp, bitscore, evalue):
        if self._is_last(sseqid):
            hit = self.hits[-1]
        else:
            hit = Hit(sseqid, slen)
            self.add(hit)
        hit.create(qstart, qend, sstart, send, mismatch, gapopen, length, pident, qcovhsp, bitscore, evalue)
        return hit

    def __str__(self):
        if self.is_empty():
            return ''
        per_hsp = []
        for hit in self.hits:
            if hit.is_empty():
                continue
            for hsp in hit.hsps:
                per_hsp.append(f'{self.qseqid}\t{self.qlen}\t{hit.sseqid}\t{hit.slen}\t{hsp}')
        return '\n'.join(per_hsp)

    def __repr__(self):
        return f'<Record for {self.qseqid}>'


class Reader:
    """Read record from BlastOutputFile object.

    Args:
        file: <BlastOutputFile>, opened BlastOutputFile object
    """
    def __init__(self, file):
        if not file.readable():
            raise IOError(f'{file.path} is not readable')
        self.file = file
        self._current = None

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            row = self.file.read()
            if row is None:
                if self._current is None:
                    raise StopIteration
                else:
                    temp = self._current
                    self._current = None
                    return temp
            if self._current is None:
                self._current = Record(row[0], row[1])
            if row[0] == self._current.qseqid:
                self._current.create(*row[2:])
            else:
                temp = self._current
                self._current = Record(row[0], row[1])
                self._current.create(*row[2:])
                return temp

    def __repr__(self):
        return f'<Reader for {self.file.path}>'


class Writer:
    """Write record to BlastOutputFile object.

    Args:
        file: <BlastOutputFile>, opened BlastOutputFile object
    """
    def __init__(self, file):
        if not file.writable():
            raise IOError(f'{file.path} is not writable')
        self.file = file

    def write(self, record):
        if record.is_empty():
            return
        for hit in record.hits:
            if hit.is_empty():
                continue
            for hsp in hit.hsps:
                self.file.write([
                    record.qseqid, record.qlen, hit.sseqid, hit.slen, hsp.qstart, hsp.qend, hsp.sstart, hsp.send,
                    hsp.mismatch, hsp.gapopen, hsp.length, hsp.pident, hsp.qcovhsp, hsp.bitscore, hsp.evalue
                ])

    def __repr__(self):
        return f'<Writer for {self.file.path}>'


class BlastOutputFile:
    """BLAST output file.

    Note: To read BLAST alignment results using the BlastOutputFile, the file the outfmt parameter must be set to
    "6 qseqid qlen sseqid slen qstart qend sstart send mismatch gapopen length pident qcovhsp bitscore evalue".

    Args:
        path: <pathlib.Path>, BLAST output file path
        mode: <str>, r: open for reading, w: open for writing
    """

    def __init__(self, path, mode='r'):
        if isinstance(path, str):
            path = Path(path)
        self.path = path
        if (mode != 'r') and (mode != 'w'):
            raise ValueError(f'invalid mode: {mode}')
        self.mode = mode
        self._fp = None

    def open(self):
        self._fp = self.path.open(self.mode)

    def close(self):
        self._fp.close()
        self._fp = None

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._fp is not None:
            self.close()

    def readable(self):
        return self._fp.readable() if self._fp else False

    def read(self):
        line = self._fp.readline()
        if not line:
            return None
        row = line.rstrip().split('\t')
        try:
            return [
                row[0], int(row[1]), row[2], int(row[3]), int(row[4]), int(row[5]), int(row[6]), int(row[7]),
                int(row[8]), int(row[9]), int(row[10]), float(row[11]), float(row[12]), float(row[13]), float(row[14])
            ]
        except (ValueError, IndexError) as exc:
            raise IOError(
                f'File {self.path} must be a tab-delimited file, and its columns must be organized in the following order: qseqid, '
                f'qlen, sseqid, slen, qstart, qend, sstart, send, mismatch, gapopen, length, pident, qcovhsp, bitscore, and evalue: {exc}'
            )

    def __iter__(self):
        if not self.readable():
            raise IOError(f'{self.path} is not readable')
        return self

    def __next__(self):
        row = self.read()
        if row is None:
            raise StopIteration
        return row

    def writable(self):
        return self._fp.writable() if self._fp else False

    def write(self, row):
        try:
            self._fp.write(
                f'{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\t{row[6]}\t{row[7]}\t'
                f'{row[8]}\t{row[9]}\t{row[10]}\t{row[11]}\t{row[12]}\t{row[13]}\t{row[14]}\n'
            )
        except IndexError:
            raise IOError(f'{row} is supposed to have 15 columns, but it only has {len(row)} columns.')

    def __repr__(self):
        return f'<BlastOutputFile {self.path}>'
