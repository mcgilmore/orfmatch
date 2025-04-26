# orfmatch

Transfer feature annotations from a reference genome to a *de novo* assembled one, where the new genome sequence is from the same or a closely related strain.

## Installation

Install using pip:

`pip install orfmatch`

or from github:

`pip install git+https://github.com/mcgilmore/orfmatch.git`

## Usage

`orfmatch [-v (Optional: outputs sequence variants as fasta and alignments)] --input <assembly.fasta> --reference <reference.gbff> --output <output.gbff>`

- Input is an assembly in \*.fasta format.
- Reference genome and output genome are in GenBank format (\*.gbff).
- Optionally, sequences which differ from the reference but are still classified as the same by pyhmmer can be output using the `-v` or `--variants` argument. A fasta file containing all varying sequences and pairwise alignments will be output to `variants.fasta` and `variants_alignment.txt` respectively.
