# orfmatch

Transfer feature annotations from a reference genome to a *de novo* assembled one, where the new genome sequence is from the same or a closely related strain.

## Installation

Install using pip:

`pip install orfmatch`

or from github:

`pip install git+https://github.com/mcgilmore/orfmatch.git`

## Usage

- Input is an assembly in \*.fasta format.
- Reference genome and output genome are in GenBank format (\*.gbff).

`orfmatch --input <assembly.fasta> --reference <reference.gbff> --output <output.gbff>`

### Optional

- `-e` / `--e-value`: E value cutoff for phmmer protein match search (default: `1e-25`).
- `-v` / `--variants`: Outputs matched sequences with differences from the reference to `variants.fasta` and alignment to `variants_alignment.txt`.
- `-c` / `--circle`: Produces a circle plot with features mapped between reference and assembly in SVG format.
- `-l` / `--line`: Produces a linear plot with features mapped between reference and assembly in SVG format.
- `-t` / `--threads`: Number of threads used for processing (default: `8`)
