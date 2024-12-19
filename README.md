# blastutils: BLAST Utilities

When using BLAST, you often need to manipulate the results, such as reading, sorting, and filtering. Biopython is a good choice, but the newer version of Biopython seems to be phasing out support for tab-delimited text files, opting for XML and XML2 formats instead. This change brings several issues:

1. XML text is not very human-readable.
2. XML files are significantly larger than tab-delimited text files. I encountered a case where the results file was 1.5TB in XML format, while the tab-delimited version was only 40GB.
3. Bio.Blast.parse has a bug that sometimes raises a CorruptedXMLError (this is critical).
4. Biopython's support for sorting, filtering, and other operations is somewhat weak.

## Installation

```shell
pip install blastutils
```

## `HSP`, `Hit`, `Record`

### `HSP`

In BLAST, an HSP (High-Scoring Segment Pair) refers to a segment of alignment between two sequences that exhibits high similarity based on certain evaluation criteria such as score, e-value, or sequence coverage. HSPs are often the most important part of a BLAST result, as they represent the strongest regions of similarity between two sequences. The HSP object has 11 attributes:

1. `qstart` – The start position of the alignment on the query sequence.
2. `qend` – The end position of the alignment on the query sequence.
3. `sstart` – The start position of the alignment on the subject sequence.
4. `send` – The end position of the alignment on the subject sequence.
5. `mismatch` – The number of mismatches in the alignment.
6. `gapopen` – The number of gaps in the alignment.
7. `length` – The length of the alignment.
8. `pident` – The percentage of identical matches.
9. `qcovhsp` – The coverage of the query sequence by this HSP.
10. `bitscore` – The bit score of the alignment.
11. `evalue` – The e-value of the alignment.

### `Hit`

In BLAST, a Hit refers to the result of an alignment between a query sequence and a subject sequence from the database. A Hit typically contains one or more HSPs. The Hit object has 3 attributes:

1. `sseqid` – The identifier of the subject sequence.
2. `slen` – The length of the subject sequence.
3. `hsps` – The HSPs associated with the hit.

### `Record`

A Record represents all the results of a single query sequence, typically containing one or more Hits. The Record object has 3 attributes:

1. `qseqid` – The identifier of the query sequence.
2. `qlen` – The length of the query sequence.
3. `hits` – The Hits associated with the record.

## Usage

### Reading and writing BLAST result file

Note: To parse BLAST alignment results using the `BlastOutputFile`, the `outfmt` parameter must be set to
`6 qseqid qlen sseqid slen qstart qend sstart send mismatch gapopen length pident qcovhsp bitscore evalue`.

```python
from blastutils import BlastOutputFile, Reader, Writer

with BlastOutputFile('example.txt') as file:
    reader = Reader(file)
    records = []
    for record in reader:
        records.append(record)

with BlastOutputFile('new-example.txt', 'w') as file:
    writer = Writer(file)
    for record in records:
        writer.write(record)
```

`BlastOutputFile` is a wrapper for BLAST result files and has two parameters:

1. `path` – The file path.
2. `mode` – The mode for opening the file; 'r' for reading (default), 'w' for writing.

You can use it with a context manager as shown above, or manage the file manually:

```python
from blastutils import BlastOutputFile, Reader

file = BlastOutputFile('example.txt')
file.open()
row = file.read()  # ['seq1', 100, 'ref1', 100, 1, 100, 1, 100, 0, 0, 100, 100.0, 100.0, 185.0, 8.22e-50]
file.close()
```

Pass the `BlastOutputFile` object to the `Reader` to read the results. The `Reader` is an iterator, and you can iterate over it using a for loop, or use next():

```python
from blastutils import BlastOutputFile, Reader

file = BlastOutputFile('example.txt')
file.open()
reader = Reader(file)
record1 = next(reader)
record2 = next(reader)
file.close()
```

### Filtering

Filtering is a very common task, for example, filtering out Hits with low similarity or low coverage. However, filtering criteria can vary widely, so `blastutils` does not provide a lot of built-in filters. Instead, it provides an API for users to implement custom filtering rules.

`blastutils` defines a `Filter` base class, and all subclasses must implement the `__call__` method, which takes a `Hit` object as input and returns `True` if the `Hit` should be kept, or `False` if it should be filtered out.

For example, to implement a filtering rule that removes Hits with an HSP's identity (pident) lower than 90:

```python
from blastutils import Filter

class MinSimilarity(Filter):
    def __init__(self, threshold):
        self.threshold = threshold

    def __call__(self, hit):
        if hit.is_empty():  # Filter out Hits with no HSPs
            return False
        return hit.hsps[0].pident >= self.threshold

record.filter(MinSimilarity(90))  # Filter out Hits with identity below 90
```

### Sorting

After obtaining BLAST results, a common task is to find the best match for each query sequence. One common approach is to sort all the Hits and select the first one. Like filtering, sorting criteria can vary, and `blastutils` does not provide a universal sorting method. Instead, it defines a Compare base class, and all subclasses must implement the `__call__` method, which takes two `Hit` objects as input:

1. Return `-1` (or any value less than `0`) if `hit1` is better than `hit2`.
2. Return `0` if `hit1` and `hit2` are equally good.
3. Return `1` (or any value greater that 0) if `hit2` is better than `hit1`.

For example, to implement a sorting rule that prefers Hits with smaller e-values and larger bit scores:

```python
from blastutils import Compare

class ByEvalueBitscore(Compare):
    def __call__(self, hit1, hit2):
        if hit2.is_empty():  # If hit2 has no HSP, hit1 is better
            return -1
        if hit1.is_empty():  # If hit1 has no HSP, hit2 is better
            return 1
        hsp1 = hit1.hsps[0]
        hsp2 = hit2.hsps[0]
        if hsp1.evalue < hsp2.evalue:
            return -1
        if hsp1.evalue > hsp2.evalue:
            return 1
        return -1 if hsp1.bitscore >= hsp2.bitscore else 1

record.best(ByEvalueBitscore())  # You can get the best Hit without sorting
record.sort(ByEvalueBitscore())  # Sort the Hits
```
