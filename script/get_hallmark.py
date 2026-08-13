#!/usr/bin/env python3
"""
get_hallmark_v3.py

Extract protein sequences for one annotated hallmark gene from a nucleotide
FASTA and GFF3 file.

Gene-name matching is CASE-SENSITIVE.

For --gene MCP, matching is performed specifically against:
    plv_annotation_db=MCP

This is intentional because MCP is curated in the supplied PLV GFF through
the plv_annotation_db field. Other gene names retain the v2 Name=/fallback
behavior.

Example:
    python get_hallmark_v3.py input.fa input.gff --gene RVE --prefix Aquinto_RVE

Output:
    Aquinto_RVE.faa

Headers are numbered and retain provenance after the first whitespace:
    RVE_001 contig=... ID=... start=... end=... strand=-

The extracted proteins do not contain a terminal stop character.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class CDS:
    seqid: str
    start: int
    end: int
    strand: str
    phase: int
    attrs: dict[str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract proteins for one hallmark gene from FASTA + GFF3."
    )
    parser.add_argument("fasta", help="Input nucleotide FASTA.")
    parser.add_argument("gff", help="Input GFF3 annotation.")
    parser.add_argument(
        "--gene",
        required=True,
        help=(
            "Gene annotation to extract. Matching is case-sensitive. "
            "For MCP, the script specifically matches plv_annotation_db=MCP."
        ),
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output prefix. Result is PREFIX.faa.",
    )
    return parser.parse_args()


def read_fasta(path: str) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current_id: str | None = None

    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n\r")
            if not line:
                continue

            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}.")
                current_id = header.split()[0]
                if current_id in sequences:
                    raise ValueError(
                        f"Duplicate FASTA sequence ID '{current_id}' at line {line_number}."
                    )
                sequences[current_id] = []
            else:
                if current_id is None:
                    raise ValueError(
                        f"FASTA sequence encountered before a header at line {line_number}."
                    )
                sequences[current_id].append("".join(line.split()).upper())

    return {key: "".join(value) for key, value in sequences.items()}


def parse_attributes(field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in field.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[unquote(key)] = unquote(value)
        else:
            attrs[unquote(item)] = ""
    return attrs


def matches_gene(attrs: dict[str, str], gene_name: str) -> bool:
    """Return True only for an exact, case-sensitive annotation match.

    MCP is special-cased because the PLV GFF uses the curated
    ``plv_annotation_db=MCP`` field for MCP classification. For MCP,
    Name= is deliberately ignored.

    All other genes keep the v2 behavior:
    Name= is the primary annotation, with plv_annotation_db= used only
    when Name= is absent.
    """
    target = gene_name.strip()
    if not target:
        return False

    if target == "MCP":
        return attrs.get("plv_annotation_db", "").strip() == "MCP"

    if "Name" in attrs:
        return attrs["Name"].strip() == target

    value = attrs.get("plv_annotation_db")
    return value is not None and value.strip() == target


def parse_gff(path: str, gene_name: str) -> list[CDS]:
    hits: list[CDS] = []

    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"GFF line {line_number} has {len(fields)} columns; expected 9."
                )

            seqid, _source, feature, start, end, _score, strand, phase, attr_text = fields
            if feature != "CDS":
                continue
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"GFF line {line_number}: unsupported strand '{strand}'."
                )

            attrs = parse_attributes(attr_text)
            if not matches_gene(attrs, gene_name):
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError as exc:
                raise ValueError(
                    f"GFF line {line_number}: start/end must be integers."
                ) from exc

            if start_i < 1 or end_i < start_i:
                raise ValueError(
                    f"GFF line {line_number}: invalid coordinates {start}-{end}."
                )

            if phase == ".":
                phase_i = 0
            else:
                try:
                    phase_i = int(phase)
                except ValueError as exc:
                    raise ValueError(
                        f"GFF line {line_number}: invalid phase '{phase}'."
                    ) from exc
                if phase_i not in {0, 1, 2}:
                    raise ValueError(
                        f"GFF line {line_number}: phase must be 0, 1, or 2."
                    )

            hits.append(CDS(seqid, start_i, end_i, strand, phase_i, attrs))

    return hits


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def translate(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    usable = len(seq) - (len(seq) % 3)
    protein = "".join(
        CODON_TABLE.get(seq[i:i + 3], "X")
        for i in range(0, usable, 3)
    )

    # Remove only a terminal stop. Internal stops are retained.
    return protein[:-1] if protein.endswith("*") else protein


def extract_single_cds(sequence: str, cds: CDS) -> str:
    if cds.end > len(sequence):
        raise ValueError(
            f"Coordinates {cds.start}-{cds.end} exceed sequence length "
            f"{len(sequence)} for '{cds.seqid}'."
        )

    # GFF3 coordinates are 1-based and inclusive.
    nucleotide = sequence[cds.start - 1:cds.end]

    if cds.strand == "-":
        nucleotide = reverse_complement(nucleotide)

    # GFF3 phase is defined from the CDS 5' end in transcriptional orientation.
    if cds.phase:
        if cds.phase >= len(nucleotide):
            raise ValueError(
                f"Phase {cds.phase} is incompatible with CDS length "
                f"{len(nucleotide)} for '{cds.seqid}:{cds.start}-{cds.end}'."
            )
        nucleotide = nucleotide[cds.phase:]

    return translate(nucleotide)


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def write_fasta(
    path: str,
    gene_name: str,
    hits: list[CDS],
    sequences: dict[str, str],
) -> int:
    count = 0
    target = gene_name.strip()

    with open(path, "w", encoding="utf-8") as out:
        for cds in sorted(hits, key=lambda x: (x.seqid, x.start, x.end)):
            sequence = sequences.get(cds.seqid)
            if sequence is None:
                raise ValueError(
                    f"GFF sequence ID '{cds.seqid}' was not found in the FASTA."
                )

            protein = extract_single_cds(sequence, cds)
            if not protein:
                continue

            count += 1
            annotation_id = (
                cds.attrs.get("ID")
                or cds.attrs.get("Name")
                or cds.attrs.get("index_name")
                or "NA"
            )

            header = (
                f"{target}_{count:03d} "
                f"contig={cds.seqid} "
                f"ID={annotation_id} "
                f"start={cds.start} "
                f"end={cds.end} "
                f"strand={cds.strand}"
            )

            out.write(f">{header}\n")
            out.write(wrap_sequence(protein) + "\n")

    return count


def main() -> int:
    args = parse_args()
    gene = args.gene.strip()

    if not gene:
        print("ERROR: --gene cannot be empty.", file=sys.stderr)
        return 1

    try:
        sequences = read_fasta(args.fasta)
        if not sequences:
            raise ValueError("No sequences were found in the FASTA.")

        hits = parse_gff(args.gff, gene)
        if not hits:
            raise ValueError(
                f"No CDS with exact case-sensitive annotation '{gene}' was found in the GFF3."
            )

        output = Path(f"{args.prefix}.faa")
        count = write_fasta(str(output), gene, hits, sequences)

        if count == 0:
            raise ValueError(
                f"Gene '{gene}' was found in the GFF3, but no protein sequence could be extracted."
            )

        print(f"Gene:    {gene}")
        print("Match:   exact, case-sensitive")
        print(f"Matched: {len(hits)}")
        print(f"Written: {count}")
        print(f"Output:  {output}")
        return 0

    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
