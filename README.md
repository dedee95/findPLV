![Python 3.10](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![GPL v3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)

# findPLV
A Robust pipeline for Polinton-like Viruses (PLV) identification in Eukaryotic Genome. It took uukaryotic genome assemblies as input, then retain PLV based on their conserved hallmark genes.

## Installation
Make sure all of these depedencies already installed in your system.
- Python 3.10+
- pyfastx
- pyhmmer
- pyrodigal
- BLAST+

Then clone this repository.
```
git clone https://github.com/dedee95/findPLV.git
```

Type -h to see full help message and to make sure all of depedencies is successfully installed.
```bash
python findPLV.py -h 

findPLV.py - Identify Polinton-Like Viruses (PLVs) in eukaryotic genome assemblies.

Usage: findPLV.py -db <directory> --prefix <prefix> <genome.fa> [OPTIONS]

Mandatory:
  -db, --db            Database directory containing:
                         PLV_hallmarks.hmm
                         MCP_nonPLV.hmm
                         Pfam-A.hmm
  --prefix             Output prefix for PLV IDs and files
  genome               Input eukaryotic genome assembly FASTA (gzip accepted)

Optionals:
  -o, --outdir         Output directory                    [default: ./Result_<YYYYMMDD>]
  -t, --threads        CPU threads                         [default: 4]
  -p, --parallel       Parallel candidate/TIR workers      [default: --threads]
  -g, --gff            Host eukaryotic GFF/GFF3 annotation. Pyrodigal ORFs
                       fully contained in host annotation intervals are removed;
                       final PLV spans overlapping host annotation are rejected.
  -e, --evalue         HMM E-value cutoff                  [default: 1e-5]
  -h, --help           Show this help and exit
```

