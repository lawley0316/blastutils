# blastutils: BLAST Utilities

## Installation

```shell
pip install blastutils
```

## Usage

### Parse and iterate over BLAST result file

Note: To parse BLAST alignment results using the `BlastOutputFile`, the `outfmt` parameter must be set to
`6 qseqid qlen sseqid slen qstart qend sstart send mismatch gapopen length pident qcovhsp bitscore evalue`.

```python
from blastutils import BlastOutputFile


with BlastOutputFile('/path/to/blast.txt') as blast_outfile:
    for record in blast_outfile:
        pass
```