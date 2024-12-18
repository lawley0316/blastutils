# blastutils: BLAST Utilities

## Installation

```shell
pip install blastutils
```

## Usage

### Read and BLAST result file

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