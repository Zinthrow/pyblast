
# PyBlast

PyBlast is a Python interface for running BLAST (Basic Local Alignment Search Tool) operations easily and efficiently. It provides simple functions to perform sequence alignments against FASTA files or prebuilt BLAST databases, with additional utilities for retrieving sequences from databases.

## Installation



To install, you can clone the repository and run pip:

```sh
git clone https://github.com/Zinthrow/pyblast.git
cd pyblast
pip install .
```

## Quick Start

First, import PyBlast:

```python
from pyblast.blastn.api import PyBlastnAPI

pybla = PyBlastnAPI()
```

### Running BLAST

You can retrieve the default BLAST parameters:

```python
default_params = pybla.wrapper.get_default_parameters()
print(default_params)  # Default parameters including 'outfmt'
```

Running BLAST against a FASTA file:

```python
outfmt = "6 qseqid sseqid pident qcovhsp length qlen mismatch gapopen qstart qend sstart send evalue bitscore"

df = pybla.blast_fasta_against_fasta(
    'segment.fasta', 'Human_RNA.fasta', outfmt=outfmt,
    num_threads=4, word_size=11
)  # Additional BLAST parameters can be passed as kwargs
```

Alternatively, if you have a prebuilt BLAST database:

```python
df = pybla.blast_fasta_against_db(
    'segment.fasta', 'human_rna_db', outfmt=outfmt,
    num_threads=4, word_size=11
)
```

### Retrieving Sequences from a Database

If you need to extract sequences that matched in a BLAST search:

```python
seqids = df.sseqid.values  # Extract matching sequence IDs
pybla.local_db.pull_fasta_from_db(seqids, 'database_hit_seqs.fasta')
```

## Installing BLAST

PyBlast provides a helper script to install BLAST from NCBI’s repository based on your operating system.

### Installation Script Usage

To install BLAST, run:

```python
from pyblast.blast_installer.install import install_blast

install_blast("2.15.0")  # Installs BLAST version 2.15.0
```

### Verifying Installation

To verify that BLAST is installed correctly:

```python
from pyblast.blast_installer.install import test_blast_installation

test_blast_installation("blastn")  # Returns True if BLAST is installed
```

### Check Database general info

```python
transcript_db = '../../software/blast_databases/blastdb/grch38_cdna'

pblast = PyBlastnAPI()
pblast.local_db.dbpath = transcript_db

pblast.local_db.get_blastdb_info()
```
```
{'Database': 'GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna',
 'Sequence Count': 145439,
 'Total Bases': 297167579,
 'Date': 'Oct 2, 2024  6:34 PM',
 'Longest Sequence': 107976,
 'BLASTDB Version': 5,
 'Volumes': '/mnt/c/Users/alars/Projects/software/blast_databases/blastdb/grch38_cdna'}
```

## File Structure

```
pyblast/
├── blastn/
│   ├── api.py          # Main API for BLAST operations
│   ├── wrapper.py      # Low-level BLAST interface
│   ├── utils.py        # Utility functions
├── blast_installer/
│   ├── install.py      # BLAST installer script
│   └── ncbi-blast-2.15.0+-x64-linux.tar.gz  # BLAST package (optional)
├── examples/
│   ├── Human_RNA.fasta
│   ├── segment.fasta
│   ├── module_testing.ipynb
├── setup.py
├── requirements.txt
└── README.md
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Feel free to open an issue or submit a pull request on GitHub.

## Author

Developed by Zinthrow.

For more information, visit [GitHub Repository](https://github.com/Zinthrow/pyblast).

