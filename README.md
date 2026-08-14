![Python 3.10](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![GPL v3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)

# findPLV
A Robust pipeline for Polinton-like Viruses (PLV) identification in eukaryotic genomes. It takes eukaryotic genome assemblies as input, then retains PLVs based on their conserved hallmark genes.

## Installation
Make sure all of these dependencies are already installed on your system.
- Python 3.10+
- pyfastx
- pyhmmer
- pyrodigal
- BLAST+

Then, clone this repository.
```
git clone https://github.com/dedee95/findPLV.git
```

`findPLV` uses 3 kinds of HMM databases: `PLV_hallmarks.hmm`, `MCP_nonPLV.hmm`, and `Pfam-A.hmm`. Make sure these 3 databases are already in the database folder. By default, `PLV_hallmarks.hmm` and `MCP_nonPLV.hmm` are already shipped with this repository. Then download the Pfam database.

```
database/
 ├── PLV_hallmarks.hmm
 ├── MCP_nonPLV.hmm
 └── Pfam-A.hmm
```

Type `-h` to see the full help message and make sure all dependencies are successfully installed.

```bash
python findPLV.py -h 

findPLV.py - Identify Polinton-Like Viruses (PLVs) in eukaryotic genome assemblies.

Usage: findPLV.py -db <directory> --prefix <prefix> <genome.fa> [OPTIONS]

Mandatory:
  -db, --db            Database directory.
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

## Usage
findPLV can takes genome assembly file in FASTA or gziped format.
```
python findPLV.py \
    -db database \
    --prefix Chlrein \
    Chlamydomonas_reinhardtii.fasta
```

A host GFF/GFF3 annotation can optionally be supplied in case you want to make sure there is no overlap between pyrodigal ORF and eukaryotic ORF.

```
python findPLV.py \
    -db database \
    --prefix Chlrein \
    -g braker.gff \
    Chlamydomonas_reinhardtii.fasta
```

A successful run produces a directory similar to:
```
Result_20260814/
├── Chlrein.summary.tsv  # Main PLV ouput summary table
├── Chlrein.func.tsv     # ORF-level annotations for every retained PLV
├── Chlrein.markerout    # A compact coordinate table containing PLVs and their genomic features
├── Chlrein.plv.fna      # FASTA file containing the nucleotide sequence of every predicted PLV
├── Chlrein.plv.pep      # FASTA file containing all proteins predicted inside the retained PLVs
├── Chlrein.plv.cds      # FASTA file containing the nucleotide CDS sequence for every ORF in the predicted PLVs
├── Chlrein.plv.gff3     # GFF3 annotation containing the complete predicted PLV
├── run.log              # Output log information
└── marker/              # Directory contains one protein FASTA for each PLV hallmark family detected in the final PLV set
    ├── Chlrein.MCP.pep
    ├── Chlrein.mCP.pep
    └── ...
``` 

## How it works
1. `findPLV` predicts protein-coding ORFs from each eligible genome contig using Pyrodigal.
2. All predicted proteins are scanned against `PLV_hallmarks.hmm` to identify PLV hallmark genes, with major capsid protein (MCP) required for PLV detection.
3. MCP-like proteins are also scanned against `MCP_nonPLV.hmm`; MCPs with stronger non-PLV evidence are excluded from PLV classification.
4. Nearby PLV hallmark genes are grouped into candidate regions. At least two distinct hallmark families, including MCP, are required.
5. Candidates containing the conserved MCP + mCP + ATPase combination are recognized as CORE3 PLV candidates.
6. Candidate boundaries are provisionally refined using hallmark positions, nearby ORFs, and local GC composition.
7. We search each candidate region for terminal inverted repeats (TIRs) using self-BLASTN. We select the best TIR pair by prioritizing retention of as many PLV hallmarks as possible.
8. If a valid TIR is found, it defines the final PLV boundary, and we search the flanking sequences for 4–20 bp target-site duplications (TSDs). If no TIR is found, only CORE3 candidates are retained.
9. Final candidates are filtered for PLV size, hallmark content, sequence quality, non-PLV MCP contamination, and optional host-GFF overlap.

The `PLV_hallmarks.hmm` file contains several conserved PLV hallmark genes collected from previously published studies: [Yutin et al. 2015](https://link.springer.com/article/10.1186/s12915-015-0207-4); [Bellas & Sommaruga 2021](https://link.springer.com/article/10.1186/s40168-020-00956-0); [Bellas et al. 2023](https://www.pnas.org/doi/10.1073/pnas.2300465120); [Bulzu et al. 2025](https://link.springer.com/article/10.1186/s40168-025-02148-0); [Bellas & Sommaruga 2026](https://www.biorxiv.org/content/biorxiv/early/2026/06/19/2026.06.19.733378.full.pdf)


```
1. pPolB: protein primed DNA polymerase
2. MCP: major capsid protein
3. mCP: minor capsid protein
4. ATPase: DNA packaging A32-like ATPase
5. YR: Tyrosine recombinase
6. RVE: retroviral integrasse
7. MTase: DNA methyltransferase
8. Tlr6FP: tetrahymena thermophila protein
9. Tet: Tet-like diocygenase
10. Endonuc: endonuclease
```