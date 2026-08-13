#!/usr/bin/env python3
"""
findPLV.py - Identify Polinton-Like Viruses (PLVs) in eukaryotic genome assemblies.

Author: Dede Kurniawan (dedekurniawan@genomics.cn)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version as _pkg_version
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import unquote

import pyfastx
import pyhmmer
import pyrodigal


# User-tunable defaults
DEFAULTS = dict(
    min_contig=5_000,
    min_plv_length=5_000,
    max_plv_length=40_000,
    seed_window=30_000,
    min_families=3,
    cluster_merge_gap=8_000,
    max_cluster_span=40_000,
    edge_gap=3_000,
    gc_window=500,
    gc_flank=10_000,
    gc_min_delta_pct=3.0,
    hmm_evalue=1e-5,
    mcp_score_margin=0.0,
    flank_for_tir=50_000,
    tir_min_insert=6_000,
    tir_max_insert=40_000,
    tir_min_len=100,
    tir_max_len=8_000,
    tir_min_id=65.0,
    tir_min_entropy=2.0,
    tir_max_kmer_frac=0.70,
    tir_max_tandem_frac=0.70,
    tir_max_tandem_period=12,
    tir_boundary_tolerance=3_000,
    tir_max_span_ratio=1.50,
    overhang_min_coding_density=0.20,
    overhang_strict_coding_density=0.35,
    overhang_max_gc_delta=5.0,
    overhang_strict_gc_delta=3.0,
    tsd_min=4,
    tsd_max=12,
    tsd_max_slide=2,
    max_n_fraction=0.05,
    dedup_min_reciprocal_overlap=0.50,
)

MCP_LABEL = "MCP"

PLV_HALLMARK_ORDER = [
    "MCP", "mCP", "ATPase", "endonuc", "MTase",
    "pPolB", "RVE", "Tet", "Tlr6FP", "YR",
]

DB_FILES = dict(
    hallmarks="PLV_hallmarks.hmm",
    mcp_nonplv="MCP_nonPLV.hmm",
    pfam="Pfam-A.hmm",
)


HELP_TEXT = """\
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
"""
USAGE_TEXT = "Usage: findPLV.py -db <DB directory> --prefix <prefix> <genome.fa> [OPTIONS]\n"


_LOG = logging.getLogger("findPLV")
OUTPUT = 25
logging.addLevelName(OUTPUT, "OUTPUT")


def _output(self, message, *args, **kwargs):
    if self.isEnabledFor(OUTPUT):
        self._log(OUTPUT, message, args, **kwargs)


logging.Logger.output = _output


def setup_logging(log_path: Optional[Path] = None) -> None:
    _LOG.handlers.clear()
    _LOG.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    sh.setLevel(logging.INFO)
    _LOG.addHandler(sh)
    if log_path is not None:
        fh = logging.FileHandler(log_path, mode="w")
        fh.setFormatter(fmt)
        fh.setLevel(logging.DEBUG)
        _LOG.addHandler(fh)


_NATKEY_RE = re.compile(r"(\d+)")


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in _NATKEY_RE.split(str(s))]


def _hallmark_sort_key(name: str):
    name = str(name)
    try:
        return (0, PLV_HALLMARK_ORDER.index(name))
    except ValueError:
        return (1, _natural_key(name))


def _fetch_seq(fa: pyfastx.Fasta, genome_path: str | Path, contig: str, start: int, end: int) -> str:
    start = max(1, int(start))
    end = max(start, int(end))
    try:
        return fa.fetch(contig, (start, end))
    except UnicodeDecodeError:
        opener = gzip.open if str(genome_path).endswith(".gz") else open
        target = contig.encode()
        seq = bytearray()
        in_target = False
        with opener(genome_path, "rb") as fh:
            for raw in fh:
                if raw.startswith(b">"):
                    name = raw[1:].strip().split(None, 1)[0]
                    if in_target:
                        break
                    in_target = name == target
                    continue
                if in_target:
                    seq.extend(raw.strip())
                    if len(seq) >= end:
                        break
        if not seq:
            raise
        return bytes(seq[start - 1:end]).decode("ascii", "replace").upper()


def gc_of_seq(seq: str) -> float:
    if not seq:
        return float("nan")
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    return 100.0 * gc / (gc + at) if (gc + at) else float("nan")


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def _overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


# =============================================================================
# Data classes
# =============================================================================
@dataclass
class Orf:
    orf_id: str
    contig: str
    start: int
    end: int
    strand: int
    protein: str
    partial: bool

    # MCP and mCP evidence and all other PLV hallmark hits.
    hallmark_hits: Dict[str, Tuple[float, float]] = field(default_factory=dict)

    family: Optional[str] = None
    family_bitscore: float = 0.0
    family_evalue: float = float("inf")
    family_source: Optional[str] = None

    plv_mcp_bitscore: float = 0.0
    plv_mcp_evalue: float = float("inf")
    nonplv_mcp_bitscore: float = 0.0
    nonplv_mcp_evalue: float = float("inf")
    best_pfam_acc: Optional[str] = None
    best_pfam_name: Optional[str] = None
    best_pfam_bitscore: float = 0.0
    best_pfam_evalue: float = float("inf")
    pfam_hits: List[Tuple[str, str, float, float]] = field(default_factory=list)


@dataclass
class TirPair:
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    tir_length: int
    insert_size: int
    tir_identity: float
    score: float
    matches: int
    total: int
    gaps: int
    tir_evalue: float = float("nan")


@dataclass
class Tsd:
    sequence_left: str
    sequence_right: str
    length: int
    mismatches: int
    identity: float
    left_shift: int
    right_shift: int


@dataclass
class Plv:
    plv_id: str
    contig: str
    contig_length: int
    start: int
    end: int
    length: int
    orfs: List[Orf]
    families_present: List[str]
    n_families: int
    mcp_best_bitscore: float
    gc_plv: float
    boundary_method: str
    tir: Optional[TirPair] = None
    tsd: Optional[Tsd] = None

    @property
    def has_tir(self) -> bool:
        return self.tir is not None

def _predict_orfs_on_contig(args: Tuple[str, str]) -> Tuple[str, List[Tuple], Optional[str]]:
    contig_name, seq = args
    try:
        gene_finder = pyrodigal.GeneFinder(meta=True)
        genes = gene_finder.find_genes(seq.encode("ascii"))
    except Exception as exc:
        return contig_name, [], f"pyrodigal failed on {contig_name}: {exc}"

    out: List[Tuple] = []
    for i, gene in enumerate(genes, start=1):
        protein = gene.translate().rstrip("*")
        if not protein:
            continue
        partial = bool(getattr(gene, "partial_begin", False) or getattr(gene, "partial_end", False))
        out.append((f"orf{i:05d}", int(gene.begin), int(gene.end), int(gene.strand), protein, partial))
    return contig_name, out, None


def predict_orfs(genome_path: Path, min_contig: int, threads: int) -> Tuple[Dict[str, Orf], Dict[str, int]]:
    fa = pyfastx.Fasta(str(genome_path), build_index=True, uppercase=True)
    work_items: List[Tuple[str, str]] = []
    contig_lengths: Dict[str, int] = {}
    n_kept = n_skipped = 0

    for rec in fa:
        clen = len(rec.seq)
        contig_lengths[rec.name] = clen
        if clen < min_contig:
            n_skipped += 1
            continue
        n_kept += 1
        work_items.append((rec.name, str(rec.seq)))

    if not work_items:
        raise RuntimeError("No contigs passed the minimum-length filter")

    executor = None
    if threads > 1 and len(work_items) > 1:
        executor = ProcessPoolExecutor(max_workers=threads)
        results_iter = executor.map(_predict_orfs_on_contig, work_items, chunksize=1)
    else:
        results_iter = (_predict_orfs_on_contig(wi) for wi in work_items)

    orfs_by_id: Dict[str, Orf] = {}
    try:
        for contig_name, records, err in results_iter:
            if err:
                _LOG.warning(f"{err}; skipping contig")
                continue
            for suffix, start, end, strand, protein, partial in records:
                oid = f"{contig_name}__{suffix}"
                orfs_by_id[oid] = Orf(
                    orf_id=oid, contig=contig_name, start=start, end=end,
                    strand=strand, protein=protein, partial=partial,
                )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    _LOG.info(
        f"ORF prediction: {len(orfs_by_id):,} ORFs on {n_kept:,} contig(s) "
        f"(>= {min_contig:,} bp; {n_skipped:,} skipped)"
    )
    if not orfs_by_id:
        raise RuntimeError("No ORFs predicted")
    return orfs_by_id, contig_lengths

_GFF_GENE_LIKE_FEATURES = {
    "gene", "mrna", "transcript", "primary_transcript",
    "lncrna", "ncrna", "rrna", "trna", "snorna", "snrna", "mirna",
    "pseudogene", "pseudogenic_transcript",
}
_GFF_PART_FEATURES = {"cds", "exon"}
_GFF_IGNORED_FEATURES = {
    "intron", "start_codon", "stop_codon", "five_prime_utr", "three_prime_utr",
    "5utr", "3utr", "utr",
}


def _parse_gff_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    attr_text = attr_text.strip()
    if not attr_text or attr_text == ".":
        return attrs
    for raw_part in attr_text.split(";"):
        part = raw_part.strip()
        if not part:
            continue
        if "=" in part:
            key, val = part.split("=", 1)
        elif " " in part:
            key, val = part.split(None, 1)
            val = val.strip().strip('"')
        else:
            continue
        key = key.strip()
        val = unquote(val.strip().strip('"'))
        if key:
            attrs[key] = val
    return attrs


def _merge_intervals(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted((min(s, e), max(s, e)) for s, e in intervals)
    merged: List[Tuple[int, int]] = []
    for s, e in ordered:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def parse_host_gff_intervals(gff_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    intervals_by_contig: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    parts_by_parent: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    n_rows = n_used = 0

    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            n_rows += 1
            contig = cols[0]
            feature = cols[2].strip().lower() if len(cols) > 2 else ""
            try:
                start, end = int(cols[3]), int(cols[4])
            except (TypeError, ValueError):
                continue
            if start <= 0 or end <= 0:
                continue
            start, end = min(start, end), max(start, end)
            if feature in _GFF_IGNORED_FEATURES:
                continue

            attrs = _parse_gff_attributes(cols[8] if len(cols) >= 9 else "")
            if feature in _GFF_GENE_LIKE_FEATURES:
                intervals_by_contig[contig].append((start, end))
                n_used += 1
                continue
            if feature in _GFF_PART_FEATURES:
                intervals_by_contig[contig].append((start, end))
                n_used += 1
                parent_text = attrs.get("Parent") or attrs.get("parent")
                if parent_text:
                    parents = [p.strip() for p in parent_text.split(",") if p.strip()]
                else:
                    fallback = attrs.get("transcript_id") or attrs.get("gene_id") or attrs.get("ID") or attrs.get("Name")
                    parents = [fallback] if fallback else []
                for parent in parents:
                    parts_by_parent[(contig, parent)].append((start, end))
                continue

            if any(tok in feature for tok in ("gene", "transcript", "mrna")):
                intervals_by_contig[contig].append((start, end))
                n_used += 1
            elif any(tok == feature or tok in feature for tok in ("cds", "exon")):
                intervals_by_contig[contig].append((start, end))
                n_used += 1

    for (contig, _parent), parts in parts_by_parent.items():
        if parts:
            intervals_by_contig[contig].append((min(s for s, _ in parts), max(e for _, e in parts)))

    merged = {contig: _merge_intervals(iv) for contig, iv in intervals_by_contig.items() if iv}
    _LOG.info(
        f"Host GFF mask: parsed {n_rows:,} feature row(s); used {n_used:,}; "
        f"collapsed to {sum(len(v) for v in merged.values()):,} interval(s) on {len(merged):,} contig(s)"
    )
    return merged


def _is_fully_contained_in_intervals(intervals: Sequence[Tuple[int, int]], start: int, end: int) -> bool:
    if not intervals:
        return False
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][0] <= start:
            lo = mid + 1
        else:
            hi = mid
    idx = lo - 1
    return idx >= 0 and intervals[idx][0] <= start and intervals[idx][1] >= end


def _span_overlaps_intervals(intervals: Sequence[Tuple[int, int]], start: int, end: int) -> bool:
    if not intervals:
        return False
    for s, e in intervals:
        if s > end:
            break
        if e >= start:
            return True
    return False


def filter_orfs_by_host_gff(orfs_by_id: Dict[str, Orf], host_intervals: Dict[str, List[Tuple[int, int]]]) -> int:
    remove_ids = [
        oid for oid, o in orfs_by_id.items()
        if _is_fully_contained_in_intervals(host_intervals.get(o.contig, ()), o.start, o.end)
    ]
    for oid in remove_ids:
        del orfs_by_id[oid]
    _LOG.info(
        f"Host GFF mask: removed {len(remove_ids):,} Pyrodigal ORF(s) fully contained in host annotations; "
        f"{len(orfs_by_id):,} ORF(s) retained"
    )
    return len(remove_ids)

def _digital_sequences(orfs: Iterable[Orf], alphabet: pyhmmer.easel.Alphabet) -> List[pyhmmer.easel.DigitalSequence]:
    seqs = []
    for o in orfs:
        try:
            t = pyhmmer.easel.TextSequence(name=o.orf_id.encode(), sequence=o.protein)
            seqs.append(t.digitize(alphabet))
        except Exception as exc:
            _LOG.warning(f"Skipping ORF {o.orf_id} during HMM digitization: {exc}")
    return seqs


def _query_meta(top_hits) -> Tuple[str, Optional[str]]:
    try:
        name_attr = top_hits.query.name
        acc_attr = top_hits.query.accession
    except AttributeError:
        name_attr = top_hits.query_name
        acc_attr = top_hits.query_accession
    name = name_attr.decode() if isinstance(name_attr, bytes) else str(name_attr)
    if acc_attr:
        acc = acc_attr.decode() if isinstance(acc_attr, bytes) else str(acc_attr)
        acc = acc.split(".")[0]
    else:
        acc = None
    return name, acc


def _hmm_profile_names(hmm_path: Path) -> List[str]:
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        hmms = list(hf)
    if not hmms:
        raise RuntimeError(f"No HMM profiles found in {hmm_path}")
    names = []
    for hmm in hmms:
        name = hmm.name
        names.append(name.decode() if isinstance(name, bytes) else str(name))
    return names


def scan_hallmarks(
    orfs_by_id: Dict[str, Orf],
    hmm_path: Path,
    evalue: float,
    threads: int,
) -> Tuple[int, int]:
    """Scan all ORFs against the comprehensive PLV hallmark HMM collection.
    """
    profile_names = _hmm_profile_names(hmm_path)
    if MCP_LABEL not in profile_names:
        raise RuntimeError(
            f"{hmm_path} does not contain an exact 'MCP' profile; "
            "MCP is mandatory for PLV detection"
        )

    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_sequences(orfs_by_id.values(), alphabet)
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        hmms = list(hf)

    n_hits = 0
    n_orfs_with_hits = set()
    for top_hits in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
        hmm_name, _ = _query_meta(top_hits)
        for hit in top_hits:
            if not hit.included:
                continue
            target = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            o = orfs_by_id.get(target)
            if o is None:
                continue
            score, ev = float(hit.score), float(hit.evalue)
            previous = o.hallmark_hits.get(hmm_name)
            if previous is None or score > previous[0]:
                o.hallmark_hits[hmm_name] = (score, ev)
            n_hits += 1
            n_orfs_with_hits.add(target)

    missing_expected = [name for name in PLV_HALLMARK_ORDER if name not in profile_names]
    if missing_expected:
        _LOG.warning(
            "PLV hallmark database is missing expected profile(s): "
            + ", ".join(missing_expected)
        )
    _LOG.info(
        f"PLV hallmark scan: {n_hits:,} included hit(s) across "
        f"{len(n_orfs_with_hits):,} ORF(s); profiles={len(profile_names):,}"
    )
    return n_hits, len(n_orfs_with_hits)


def scan_nonplv_mcp(
    orfs_by_id: Dict[str, Orf],
    hmm_path: Path,
    evalue: float,
    threads: int,
) -> int:
    """Scan ORFs against the NCLDV/adenovirus MCP competitor collection."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_sequences(orfs_by_id.values(), alphabet)
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        hmms = list(hf)
    if not hmms:
        raise RuntimeError(f"No HMM profiles found in {hmm_path}")

    n_hits = 0
    for top_hits in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
        for hit in top_hits:
            if not hit.included:
                continue
            target = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            o = orfs_by_id.get(target)
            if o is None:
                continue
            score, ev = float(hit.score), float(hit.evalue)
            if score > o.nonplv_mcp_bitscore:
                o.nonplv_mcp_bitscore = score
                o.nonplv_mcp_evalue = ev
            n_hits += 1

    _LOG.info(
        f"{hmm_path.name}: {n_hits:,} included competitor MCP hit(s) "
        f"at E-value <= {evalue:g}"
    )
    return n_hits


def _best_hallmark_except(o: Orf, excluded: Sequence[str]) -> Tuple[Optional[str], float, float]:
    excluded_set = set(excluded)
    candidates = [
        (name, score, ev)
        for name, (score, ev) in o.hallmark_hits.items()
        if name not in excluded_set
    ]
    if not candidates:
        return None, 0.0, float("inf")
    candidates.sort(key=lambda x: (x[1], -x[2]), reverse=True)
    return candidates[0]


def classify_hallmarks(
    orfs_by_id: Dict[str, Orf],
    mcp_score_margin: float,
) -> Tuple[int, int, int]:
    """Assign one primary PLV hallmark per ORF after MCP competition.
    """
    n_mcp = n_ambiguous = n_marked = 0
    for o in orfs_by_id.values():
        o.family = None
        o.family_bitscore = 0.0
        o.family_evalue = float("inf")
        o.family_source = None
        o.plv_mcp_bitscore = 0.0
        o.plv_mcp_evalue = float("inf")
        mcp = o.hallmark_hits.get(MCP_LABEL)
        if mcp is not None:
            o.plv_mcp_bitscore, o.plv_mcp_evalue = mcp
            if o.nonplv_mcp_bitscore <= 0.0 or mcp[0] > o.nonplv_mcp_bitscore + mcp_score_margin:
                o.family = MCP_LABEL
                o.family_bitscore = mcp[0]
                o.family_evalue = mcp[1]
                o.family_source = "PLV_hallmark"
                n_mcp += 1
            else:
                n_ambiguous += 1

        if o.family is None:
            name, score, ev = _best_hallmark_except(o, (MCP_LABEL,))
            if name is not None:
                o.family = name
                o.family_bitscore = score
                o.family_evalue = ev
                o.family_source = "PLV_hallmark"
                n_marked += 1
        elif o.family == MCP_LABEL:
            n_marked += 1

    _LOG.info(
        f"Competitive PLV hallmark classification: MCP={n_mcp:,}; "
        f"competitor-dominated MCP={n_ambiguous:,}; other PLV hallmark ORFs={n_marked:,}"
    )
    return n_mcp, n_ambiguous, n_marked


def scan_pfam(
    orfs_to_scan: List[Orf],
    pfam_hmm_path: Path,
    evalue: float,
    threads: int,
) -> None:
    """Annotate candidate-region proteins with their best Pfam hit.
    """
    if not orfs_to_scan:
        return
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_sequences(orfs_to_scan, alphabet)
    by_id = {o.orf_id: o for o in orfs_to_scan}
    n_hits = 0

    with pyhmmer.plan7.HMMFile(str(pfam_hmm_path)) as hf:
        for top_hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, E=evalue):
            hmm_name, hmm_acc = _query_meta(top_hits)
            acc_rec = hmm_acc or hmm_name
            for hit in top_hits:
                if not hit.included:
                    continue
                target = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                o = by_id.get(target)
                if o is None:
                    continue
                score, ev = float(hit.score), float(hit.evalue)
                o.pfam_hits.append((acc_rec, hmm_name, score, ev))
                if score > o.best_pfam_bitscore:
                    o.best_pfam_acc = acc_rec
                    o.best_pfam_name = hmm_name
                    o.best_pfam_bitscore = score
                    o.best_pfam_evalue = ev
                n_hits += 1

    _LOG.info(f"Pfam-A scan: {n_hits:,} included hit(s) for candidate-region proteins")

def _is_plv_hallmark(o: Orf) -> bool:
    return o.family_source == "PLV_hallmark" and o.family is not None


def find_seed_regions(
    orfs_by_id: Dict[str, Orf],
    window_size: int,
    min_families: int,
    cluster_merge_gap: int,
    max_cluster_span: int,
) -> List[dict]:
    by_contig: Dict[str, List[Orf]] = defaultdict(list)
    for o in orfs_by_id.values():
        by_contig[o.contig].append(o)
    for orfs in by_contig.values():
        orfs.sort(key=lambda x: x.start)

    half = window_size // 2
    raw: List[dict] = []
    for contig, orfs in by_contig.items():
        markers = [o for o in orfs if _is_plv_hallmark(o)]
        if not markers:
            continue
        for anchor in markers:
            mid = (anchor.start + anchor.end) // 2
            ws, we = mid - half, mid + half
            members = [o for o in markers if o.start <= we and o.end >= ws]
            fams = {o.family for o in members if o.family}
            if len(fams) < min_families or MCP_LABEL not in fams:
                continue
            raw.append(dict(
                contig=contig,
                cluster_start=min(o.start for o in members),
                cluster_end=max(o.end for o in members),
                orf_ids=sorted({o.orf_id for o in members}),
                families=sorted(fams, key=_hallmark_sort_key),
            ))

    raw.sort(key=lambda c: (c["contig"], c["cluster_start"], c["cluster_end"]))
    pass1: List[dict] = []
    for c in raw:
        if pass1 and pass1[-1]["contig"] == c["contig"] and c["cluster_start"] <= pass1[-1]["cluster_end"]:
            m = pass1[-1]
            m["cluster_end"] = max(m["cluster_end"], c["cluster_end"])
            m["orf_ids"] = sorted(set(m["orf_ids"]) | set(c["orf_ids"]))
            m["families"] = sorted(set(m["families"]) | set(c["families"]), key=_hallmark_sort_key)
        else:
            pass1.append(dict(c))

    merged: List[dict] = []
    for c in pass1:
        if merged and merged[-1]["contig"] == c["contig"]:
            gap = c["cluster_start"] - merged[-1]["cluster_end"] - 1
            span = c["cluster_end"] - merged[-1]["cluster_start"] + 1
            if gap <= cluster_merge_gap and span <= max_cluster_span:
                m = merged[-1]
                m["cluster_end"] = max(m["cluster_end"], c["cluster_end"])
                m["orf_ids"] = sorted(set(m["orf_ids"]) | set(c["orf_ids"]))
                m["families"] = sorted(set(m["families"]) | set(c["families"]), key=_hallmark_sort_key)
                continue
        merged.append(dict(c))

    merged = [c for c in merged if MCP_LABEL in c["families"] and len(c["families"]) >= min_families]
    _LOG.info(
        f"Seeding: {len(merged):,} candidate PLV region(s) "
        f"[window={window_size:,} bp; min_hallmark_families={min_families}; "
        f"merge_gap={cluster_merge_gap:,} bp; max_span={max_cluster_span:,} bp]"
    )
    return merged

def _coding_density(orfs: Sequence[Orf], start: int, end: int) -> float:
    span = end - start + 1
    if span <= 0:
        return 0.0
    iv = []
    for o in orfs:
        s, e = max(start, o.start), min(end, o.end)
        if e >= s:
            iv.append((s, e))
    merged = _merge_intervals(iv)
    coding = sum(e - s + 1 for s, e in merged)
    return coding / span


def marker_gc_boundary(
    contig_orfs: List[Orf],
    cluster_start: int,
    cluster_end: int,
    contig_length: int,
    fa: pyfastx.Fasta,
    genome_path: str,
    cfg: dict,
) -> Tuple[int, int, List[str]]:
    notes: List[str] = []
    cluster_orfs = [o for o in contig_orfs if o.start <= cluster_end and o.end >= cluster_start]
    anchors = [o for o in cluster_orfs if _is_plv_hallmark(o)]
    if not anchors:
        return cluster_start, cluster_end, ["No anchors available; using seed-cluster boundary"]

    anchors.sort(key=lambda o: o.start)
    anchor_min = anchors[0].start
    anchor_max = max(o.end for o in anchors)
    start, end = anchor_min, anchor_max
    notes.append(f"Anchor boundary {start:,}-{end:,}")

    all_sorted = sorted(contig_orfs, key=lambda o: o.start)
    left_candidates = [o for o in all_sorted if o.end < start and o.end >= cluster_start - cfg["edge_gap"]]
    for o in reversed(left_candidates):
        gap = start - o.end - 1
        if gap > cfg["edge_gap"]:
            break
        start = o.start

    right_candidates = [o for o in all_sorted if o.start > end and o.start <= cluster_end + cfg["edge_gap"]]
    for o in right_candidates:
        gap = o.start - end - 1
        if gap > cfg["edge_gap"]:
            break
        end = o.end

    # Hard length cap while preserving every anchor.
    if end - start + 1 > cfg["max_plv_length"]:
        slack = cfg["max_plv_length"] - (anchor_max - anchor_min + 1)
        if slack <= 0:
            start, end = anchor_min, anchor_max
            notes.append("Anchor span exceeds max_plv_length; preserved anchors only")
        else:
            left_space = anchor_min - start
            right_space = end - anchor_max
            left_take = min(left_space, slack // 2)
            right_take = min(right_space, slack - left_take)
            leftover = slack - left_take - right_take
            if leftover > 0:
                extra_left = min(left_space - left_take, leftover)
                left_take += extra_left
                leftover -= extra_left
            if leftover > 0:
                right_take += min(right_space - right_take, leftover)
            start, end = anchor_min - left_take, anchor_max + right_take
            notes.append(f"Length cap applied: {start:,}-{end:,}")

    # Conservative GC refinement. It only shrinks non-anchor flanks.
    core_gc = gc_of_seq(_fetch_seq(fa, genome_path, contig_orfs[0].contig, anchor_min, anchor_max))
    if not math.isnan(core_gc):
        flank = cfg["gc_flank"]
        win = cfg["gc_window"]
        left_host_gc = gc_of_seq(_fetch_seq(fa, genome_path, contig_orfs[0].contig, max(1, start - flank), max(1, start - 1))) if start > 1 else float("nan")
        right_host_gc = gc_of_seq(_fetch_seq(fa, genome_path, contig_orfs[0].contig, min(contig_length, end + 1), min(contig_length, end + flank))) if end < contig_length else float("nan")

        if not math.isnan(left_host_gc) and abs(core_gc - left_host_gc) >= cfg["gc_min_delta_pct"]:
            pos = start
            while pos + win - 1 < anchor_min:
                wgc = gc_of_seq(_fetch_seq(fa, genome_path, contig_orfs[0].contig, pos, pos + win - 1))
                if math.isnan(wgc) or abs(wgc - core_gc) <= abs(wgc - left_host_gc):
                    break
                pos += win
            if pos > start:
                notes.append(f"GC refinement trimmed left edge {start:,}->{pos:,}")
                start = pos

        if not math.isnan(right_host_gc) and abs(core_gc - right_host_gc) >= cfg["gc_min_delta_pct"]:
            pos = end
            while pos - win + 1 > anchor_max:
                wgc = gc_of_seq(_fetch_seq(fa, genome_path, contig_orfs[0].contig, pos - win + 1, pos))
                if math.isnan(wgc) or abs(wgc - core_gc) <= abs(wgc - right_host_gc):
                    break
                pos -= win
            if pos < end:
                notes.append(f"GC refinement trimmed right edge {end:,}->{pos:,}")
                end = pos

    return max(1, start), min(contig_length, end), notes

_BLASTN_OUTFMT = "6 qstart qend sstart send length nident pident gaps evalue bitscore"


def parse_blastn_tabular(tab_path: Path) -> List[TirPair]:
    if not tab_path.exists() or not tab_path.read_text().strip():
        return []
    pairs: List[TirPair] = []
    for line in tab_path.read_text().splitlines():
        fields = line.strip().split("\t")
        if len(fields) < 10:
            continue
        try:
            qstart, qend, sstart, send = map(int, fields[:4])
            aln_len = int(fields[4])
            nident = int(fields[5])
            pident = float(fields[6])
            gaps = int(fields[7])
            evalue = float(fields[8])
            bitscore = float(fields[9])
        except ValueError:
            continue
        ls, le = min(qstart, qend), max(qstart, qend)
        rs, re = min(sstart, send), max(sstart, send)
        if le >= rs:
            continue
        tir_len = le - ls + 1
        insert = rs - le - 1
        pairs.append(TirPair(
            left_start=ls, left_end=le, right_start=rs, right_end=re,
            tir_length=tir_len, insert_size=insert, tir_identity=pident,
            score=bitscore, matches=nident, total=aln_len, gaps=gaps,
            tir_evalue=evalue,
        ))
    return pairs


def run_blastn_self(region_fa: Path, tab_out: Path) -> None:
    cmd = [
        "blastn", "-query", str(region_fa), "-subject", str(region_fa),
        "-strand", "minus", "-task", "blastn", "-word_size", "7",
        "-reward", "1", "-penalty", "-1", "-gapopen", "2", "-gapextend", "1",
        "-evalue", "10.0", "-dust", "no", "-soft_masking", "false",
        "-max_target_seqs", "10000", "-num_threads", "1",
        "-outfmt", _BLASTN_OUTFMT, "-out", str(tab_out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr)


def _dinucleotide_entropy(seq: str) -> float:
    seq = seq.upper()
    counts = Counter(seq[i:i + 2] for i in range(len(seq) - 1) if "N" not in seq[i:i + 2])
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * math.log2(c / total) for c in counts.values())


def _max_kmer_frac(seq: str, k: int) -> float:
    seq = seq.upper()
    if len(seq) < k:
        return 0.0
    best = 0.0
    for phase in range(k):
        kmers = [
            seq[i:i + k] for i in range(phase, len(seq) - k + 1, k)
            if "N" not in seq[i:i + k]
        ]
        if not kmers:
            continue
        counts = Counter(kmers)
        best = max(best, max(counts.values()) / len(kmers))
    return best


def _max_tandem_period_fraction(seq: str, max_period: int) -> float:
    seq = seq.upper()
    best = 0.0
    for p in range(1, min(max_period, len(seq) - 1) + 1):
        valid = [(a, b) for a, b in zip(seq[:-p], seq[p:]) if a != "N" and b != "N"]
        if not valid:
            continue
        frac = sum(a == b for a, b in valid) / len(valid)
        best = max(best, frac)
    return best


def _tir_is_low_complexity(seq: str, cfg: dict) -> bool:
    if not seq or _dinucleotide_entropy(seq) < cfg["tir_min_entropy"]:
        return True
    if any(_max_kmer_frac(seq, k) > cfg["tir_max_kmer_frac"] for k in (1, 2, 3, 4)):
        return True
    return _max_tandem_period_fraction(seq, cfg["tir_max_tandem_period"]) > cfg["tir_max_tandem_frac"]


def _count_bracketed(tir: TirPair, intervals: Sequence[Tuple[int, int]]) -> int:
    return sum(1 for s, e in intervals if tir.left_start <= s and e <= tir.right_end)


def select_best_tir(
    pairs: List[TirPair],
    region_offset: int,
    marker_intervals: List[Tuple[int, int]],
    region_seq: str,
    cfg: dict,
) -> Tuple[Optional[TirPair], dict]:
    diag = dict(raw=len(pairs), pass_size=0, pass_len=0, pass_id=0, pass_complexity=0, pass_bracket=0, best_near_miss="")
    valid: List[Tuple[int, TirPair]] = []
    near_score = -1.0
    required = len(marker_intervals)

    for t in pairs:
        gls, gle = t.left_start + region_offset - 1, t.left_end + region_offset - 1
        grs, gre = t.right_start + region_offset - 1, t.right_end + region_offset - 1
        abs_tir = TirPair(
            left_start=gls, left_end=gle, right_start=grs, right_end=gre,
            tir_length=t.tir_length, insert_size=t.insert_size, tir_identity=t.tir_identity,
            score=t.score, matches=t.matches, total=t.total, gaps=t.gaps, tir_evalue=t.tir_evalue,
        )
        ok_size = cfg["tir_min_insert"] <= t.insert_size <= cfg["tir_max_insert"]
        ok_len = cfg["tir_min_len"] <= t.tir_length <= cfg["tir_max_len"]
        ok_id = t.tir_identity >= cfg["tir_min_id"]
        left_seq = region_seq[t.left_start - 1:t.left_end]
        right_seq = region_seq[t.right_start - 1:t.right_end]
        ok_complex = not (_tir_is_low_complexity(left_seq, cfg) or _tir_is_low_complexity(right_seq, cfg))
        bracketed = _count_bracketed(abs_tir, marker_intervals)
        ok_bracket = bracketed >= required if required else True

        if ok_size:
            diag["pass_size"] += 1
        if ok_size and ok_len:
            diag["pass_len"] += 1
        if ok_size and ok_len and ok_id:
            diag["pass_id"] += 1
        if ok_size and ok_len and ok_id and ok_complex:
            diag["pass_complexity"] += 1
        if ok_size and ok_len and ok_id and ok_complex and ok_bracket:
            diag["pass_bracket"] += 1

        score = 1e6 * sum((ok_size, ok_len, ok_id, ok_complex, ok_bracket)) + bracketed * 1e4 + t.tir_identity * t.tir_length
        if score > near_score:
            near_score = score
            diag["best_near_miss"] = (
                f"insert={t.insert_size:,},len={t.tir_length},id={t.tir_identity:.1f}%,"
                f"bracket={bracketed}/{required},complex={'pass' if ok_complex else 'fail'}"
            )
        if ok_size and ok_len and ok_id and ok_complex and ok_bracket:
            valid.append((bracketed, abs_tir))

    if not valid:
        return None, diag
    valid.sort(key=lambda x: (x[0], x[1].tir_identity, x[1].tir_length, x[1].score), reverse=True)
    return valid[0][1], diag

def _overhang_supported(
    contig: str,
    start: int,
    end: int,
    core_gc: float,
    contig_orfs: Sequence[Orf],
    host_intervals: Sequence[Tuple[int, int]],
    fa: pyfastx.Fasta,
    genome_path: str,
    cfg: dict,
    strict: bool,
) -> bool:
    if end < start:
        return True
    if _span_overlaps_intervals(host_intervals, start, end):
        return False
    orfs = [o for o in contig_orfs if o.start <= end and o.end >= start]
    if not orfs:
        return False
    marker_present = any(_is_plv_hallmark(o) for o in orfs)
    density = _coding_density(orfs, start, end)
    min_density = cfg["overhang_strict_coding_density"] if strict else cfg["overhang_min_coding_density"]
    if density < min_density and not marker_present:
        return False
    over_gc = gc_of_seq(_fetch_seq(fa, genome_path, contig, start, end))
    max_delta = cfg["overhang_strict_gc_delta"] if strict else cfg["overhang_max_gc_delta"]
    if not math.isnan(core_gc) and not math.isnan(over_gc) and abs(core_gc - over_gc) > max_delta and not marker_present:
        return False
    return True


def arbitrate_final_boundary(
    contig: str,
    bio_start: int,
    bio_end: int,
    best_tir: Optional[TirPair],
    contig_orfs: Sequence[Orf],
    host_intervals: Sequence[Tuple[int, int]],
    fa: pyfastx.Fasta,
    genome_path: str,
    cfg: dict,
) -> Tuple[int, int, str, Optional[TirPair], str]:
    if best_tir is None:
        return bio_start, bio_end, "marker_gc_boundary", None, "no_tir"

    tol = cfg["tir_boundary_tolerance"]
    left_diff = abs(best_tir.left_start - bio_start)
    right_diff = abs(best_tir.right_end - bio_end)
    if left_diff <= tol and right_diff <= tol:
        return best_tir.left_start, best_tir.right_end, "TIR", best_tir, "boundary_compatible"

    if best_tir.left_start >= bio_start and best_tir.right_end <= bio_end:
        return bio_start, bio_end, "marker_gc_boundary", None, "internal_tir_rejected"

    # Mixed inside/outside geometry is treated conservatively.
    brackets_bio = best_tir.left_start <= bio_start and best_tir.right_end >= bio_end
    if not brackets_bio:
        return bio_start, bio_end, "marker_gc_boundary", None, "discordant_tir_rejected"

    bio_len = max(1, bio_end - bio_start + 1)
    tir_len = max(1, best_tir.right_end - best_tir.left_start + 1)
    strict = (tir_len / bio_len) > cfg["tir_max_span_ratio"]
    core_gc = gc_of_seq(_fetch_seq(fa, genome_path, contig, bio_start, bio_end))
    left_ok = _overhang_supported(
        contig, best_tir.left_start, bio_start - 1, core_gc, contig_orfs,
        host_intervals, fa, genome_path, cfg, strict,
    )
    right_ok = _overhang_supported(
        contig, bio_end + 1, best_tir.right_end, core_gc, contig_orfs,
        host_intervals, fa, genome_path, cfg, strict,
    )
    if left_ok and right_ok:
        return best_tir.left_start, best_tir.right_end, "TIR", best_tir, "outside_supported"
    return bio_start, bio_end, "marker_gc_boundary", None, "overextended_tir_rejected"


def find_tsd(left_flank: str, right_flank: str, k_min: int, k_max: int, max_slide: int) -> Optional[Tsd]:
    left, right = left_flank.upper(), right_flank.upper()
    best: Optional[Tsd] = None
    for k in range(k_max, k_min - 1, -1):
        if k > len(left) or k > len(right):
            continue
        max_mm = 0 if k <= 5 else (1 if k <= 8 else 2)
        for sl in range(max_slide + 1):
            if k + sl > len(left):
                break
            for sr in range(max_slide + 1):
                if k + sr > len(right):
                    break
                lk = left[len(left) - k - sl:len(left) - sl] if sl else left[-k:]
                rk = right[sr:sr + k]
                if "N" in lk or "N" in rk:
                    continue
                mm = sum(a != b for a, b in zip(lk, rk))
                if mm <= max_mm:
                    cand = Tsd(lk, rk, k, mm, 100.0 * (k - mm) / k, sl, sr)
                    if best is None or (cand.length, cand.identity) > (best.length, best.identity):
                        best = cand
        if best is not None and best.length == k:
            return best
    return best

def _process_cluster(task: dict) -> dict:
    ci = task["cluster_index"]
    contig = task["contig"]
    cstart, cend = task["cluster_start"], task["cluster_end"]
    clen = task["contig_length"]
    cfg = task["cfg"]
    genome_path = task["genome_path"]
    contig_orfs: List[Orf] = task["contig_orfs"]
    host_intervals: List[Tuple[int, int]] = task.get("host_intervals", [])
    fa = pyfastx.Fasta(str(genome_path), build_index=True, uppercase=True)
    logs: List[Tuple[str, str]] = []

    bio_start, bio_end, notes = marker_gc_boundary(
        contig_orfs, cstart, cend, clen, fa, genome_path, cfg,
    )
    for note in notes:
        logs.append(("info", f"Cluster {ci} {contig}:{cstart:,}-{cend:,}: {note}"))

    marker_intervals = [
        (o.start, o.end) for o in contig_orfs
        if _is_plv_hallmark(o) and o.start >= bio_start and o.end <= bio_end
    ]
    rstart = max(1, bio_start - cfg["flank_for_tir"])
    rend = min(clen, bio_end + cfg["flank_for_tir"])
    region_seq = _fetch_seq(fa, genome_path, contig, rstart, rend)

    best_tir: Optional[TirPair] = None
    with tempfile.TemporaryDirectory(prefix="findPLV_v3_") as td:
        tdir = Path(td)
        region_fa = tdir / "region.fa"
        blast_tab = tdir / "tir.tsv"
        with open(region_fa, "w") as fh:
            fh.write(f">{contig}_{rstart}_{rend}\n")
            for i in range(0, len(region_seq), 80):
                fh.write(region_seq[i:i + 80] + "\n")
        try:
            run_blastn_self(region_fa, blast_tab)
            pairs = parse_blastn_tabular(blast_tab)
            best_tir, diag = select_best_tir(pairs, rstart, marker_intervals, region_seq, cfg)
            if best_tir is None:
                logs.append((
                    "info",
                    f"Cluster {ci}: no accepted TIR [raw={diag['raw']},size={diag['pass_size']},"
                    f"len={diag['pass_len']},id={diag['pass_id']},complex={diag['pass_complexity']},"
                    f"bracket={diag['pass_bracket']}; near={diag['best_near_miss'] or 'none'}]",
                ))
        except (OSError, subprocess.CalledProcessError) as exc:
            logs.append(("warning", f"Cluster {ci}: BLASTN TIR search failed ({exc}); continuing TIR-less"))

    final_start, final_end, boundary_method, tir_for_output, tir_status = arbitrate_final_boundary(
        contig, bio_start, bio_end, best_tir, contig_orfs, host_intervals,
        fa, genome_path, cfg,
    )
    logs.append((
        "info",
        f"Cluster {ci}: boundary={boundary_method}; TIR_status={tir_status}; "
        f"biological={bio_start:,}-{bio_end:,}; final={final_start:,}-{final_end:,}",
    ))

    length = final_end - final_start + 1
    if length < cfg["min_plv_length"] or length > cfg["max_plv_length"]:
        return dict(status="skip", cluster_index=ci, log_msgs=logs,
                    message=f"Cluster {ci}: final length {length:,} bp outside {cfg['min_plv_length']:,}-{cfg['max_plv_length']:,} bp")

    # When GFF is supplied, a final PLV must not overlap any host annotation span.
    if host_intervals and _span_overlaps_intervals(host_intervals, final_start, final_end):
        return dict(status="skip", cluster_index=ci, log_msgs=logs,
                    message=f"Cluster {ci}: final PLV span overlaps host GFF annotation; discarded")

    plv_orfs = sorted(
        [o for o in contig_orfs if o.start >= final_start and o.end <= final_end],
        key=lambda o: o.start,
    )
    fams = sorted({o.family for o in plv_orfs if _is_plv_hallmark(o)}, key=_natural_key)
    if MCP_LABEL not in fams:
        return dict(status="skip", cluster_index=ci, log_msgs=logs,
                    message=f"Cluster {ci}: no competitively classified PLV MCP inside final boundary")
    if len(fams) < cfg["min_families"]:
        return dict(status="skip", cluster_index=ci, log_msgs=logs,
                    message=f"Cluster {ci}: only {len(fams)} PLV marker families remain inside final boundary")

    seq = _fetch_seq(fa, genome_path, contig, final_start, final_end)
    n_fraction = seq.upper().count("N") / len(seq) if seq else 1.0
    if n_fraction > cfg["max_n_fraction"]:
        return dict(status="skip", cluster_index=ci, log_msgs=logs,
                    message=f"Cluster {ci}: N fraction {n_fraction:.2%} > {cfg['max_n_fraction']:.0%}")

    tsd = None
    if tir_for_output is not None:
        left_flank = _fetch_seq(fa, genome_path, contig, max(1, tir_for_output.left_start - 50), tir_for_output.left_start - 1) if tir_for_output.left_start > 1 else ""
        right_flank = _fetch_seq(fa, genome_path, contig, tir_for_output.right_end + 1, min(clen, tir_for_output.right_end + 50)) if tir_for_output.right_end < clen else ""
        tsd = find_tsd(left_flank, right_flank, cfg["tsd_min"], cfg["tsd_max"], cfg["tsd_max_slide"])

    mcp_best = max((o.plv_mcp_bitscore for o in plv_orfs if o.family == MCP_LABEL), default=0.0)
    plv = Plv(
        plv_id="TBD", contig=contig, contig_length=clen,
        start=final_start, end=final_end, length=length, orfs=plv_orfs,
        families_present=fams, n_families=len(fams), mcp_best_bitscore=mcp_best,
        gc_plv=gc_of_seq(seq), boundary_method=boundary_method,
        tir=tir_for_output, tsd=tsd,
    )
    return dict(status="ok", cluster_index=ci, log_msgs=logs, plv=plv)


# =============================================================================
# Stage 8: final overlap resolution
# =============================================================================
def deduplicate_plvs(plvs: List[Plv], min_reciprocal_overlap: float) -> List[Plv]:
    if len(plvs) <= 1:
        return plvs
    ranked = sorted(
        plvs,
        key=lambda p: (p.has_tir, p.n_families, p.mcp_best_bitscore, -p.length),
        reverse=True,
    )
    kept: List[Plv] = []
    for cand in ranked:
        redundant = False
        for good in kept:
            if cand.contig != good.contig:
                continue
            overlap = _overlap_len(cand.start, cand.end, good.start, good.end)
            shorter = min(cand.length, good.length)
            if shorter and overlap / shorter >= min_reciprocal_overlap:
                _LOG.info(
                    f"Deduplication: dropped {cand.contig}:{cand.start:,}-{cand.end:,}; "
                    f"{100 * overlap / shorter:.0f}% overlap with stronger call {good.start:,}-{good.end:,}"
                )
                redundant = True
                break
        if not redundant:
            kept.append(cand)
    return kept


# =============================================================================
# Output writers
# =============================================================================
def _tsd_status(tsd: Optional[Tsd]) -> str:
    if tsd is None:
        return "NODETECT"
    return "PERFECT" if tsd.mismatches == 0 else "IMPERFECT"


def _tir_fields(tir: Optional[TirPair]) -> dict:
    if tir is None:
        return dict(
            tir_length="NA", tir_score="NA",
            tir_identity_pct="NA", tir_gaps="NA",
        )
    return dict(
        tir_length=tir.tir_length,
        tir_score=f"{tir.score:.1f}",
        tir_identity_pct=f"{tir.tir_identity:.2f}",
        tir_gaps=tir.gaps,
    )


def _tsd_fields(tsd: Optional[Tsd]) -> dict:
    if tsd is None:
        return dict(
            tsd_len="NA", tsd_left="NA", tsd_right="NA",
            tsd_mismatch="NA", tsd_conservation="NODETECT",
        )
    return dict(
        tsd_len=tsd.length,
        tsd_left=tsd.sequence_left,
        tsd_right=tsd.sequence_right,
        tsd_mismatch=tsd.mismatches,
        tsd_conservation=_tsd_status(tsd),
    )


def calculate_confidence(p: Plv) -> Tuple[int, str]:
    """Rank retained PLVs using orthogonal evidence, not raw HMM scores.
    """
    if MCP_LABEL not in p.families_present:
        return 0, "LOW"

    score = 50
    supporting = len([f for f in p.families_present if f != MCP_LABEL])
    score += min(supporting, 3) * 10
    if p.has_tir:
        score += 15
    if p.tsd is not None:
        score += 5

    score = min(100, score)
    if score >= 85:
        label = "HIGH"
    elif score >= 70:
        label = "MEDIUM"
    else:
        label = "LOW"
    return score, label


def write_summary_tsv(plvs: List[Plv], path: Path) -> None:
    columns = [
        "contig_id", "plv_name", "start", "end", "plv_length", "gc",
        "total_cds", "n_hallmark_families", "hallmarks", "mcp_bitscore",
        "confidence_score", "confidence", "has_tir",
        "tir_length", "tir_score", "tir_identity_pct", "tir_gaps",
        "tsd_len", "tsd_left", "tsd_right", "tsd_mismatch",
        "tsd_conservation", "boundary_method",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for p in sorted(plvs, key=lambda x: _natural_key(x.plv_id)):
            confidence_score, confidence = calculate_confidence(p)
            row = dict(
                contig_id=p.contig,
                plv_name=p.plv_id,
                start=p.start,
                end=p.end,
                plv_length=p.length,
                gc="NA" if math.isnan(p.gc_plv) else f"{p.gc_plv:.2f}",
                total_cds=len(p.orfs),
                n_hallmark_families=p.n_families,
                hallmarks=",".join(p.families_present),
                mcp_bitscore=f"{p.mcp_best_bitscore:.1f}",
                confidence_score=confidence_score,
                confidence=confidence,
                has_tir="yes" if p.has_tir else "no",
                boundary_method=p.boundary_method,
            )
            row.update(_tir_fields(p.tir))
            row.update(_tsd_fields(p.tsd))
            writer.writerow(row)


def _annotation_name(o: Orf) -> str:
    """Return the GFF/inspection name: hallmark first, then Pfam, then hypothetical."""
    if _is_plv_hallmark(o):
        return o.family or "hypothetical"
    if o.best_pfam_name:
        return o.best_pfam_name
    return "hypothetical"


def write_func_tsv(plvs: List[Plv], path: Path) -> None:
    columns = [
        "plv_name", "orf", "start", "end", "strand", "name", "hallmark",
        "hallmark_bitscore", "hallmark_evalue", "plv_mcp_bitscore",
        "nonplv_mcp_bitscore", "pfam_id", "pfam_name", "pfam_bitscore", "pfam_evalue",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for p in sorted(plvs, key=lambda x: _natural_key(x.plv_id)):
            for i, o in enumerate(p.orfs, start=1):
                writer.writerow(dict(
                    plv_name=p.plv_id,
                    orf=f"orf{i:05d}",
                    start=o.start,
                    end=o.end,
                    strand="+" if o.strand >= 0 else "-",
                    name=_annotation_name(o),
                    hallmark=o.family if _is_plv_hallmark(o) else "NA",
                    hallmark_bitscore=f"{o.family_bitscore:.1f}" if _is_plv_hallmark(o) and o.family_bitscore else "NA",
                    hallmark_evalue=f"{o.family_evalue:.2e}" if _is_plv_hallmark(o) and math.isfinite(o.family_evalue) else "NA",
                    plv_mcp_bitscore=f"{o.plv_mcp_bitscore:.1f}" if o.plv_mcp_bitscore else "NA",
                    nonplv_mcp_bitscore=f"{o.nonplv_mcp_bitscore:.1f}" if o.nonplv_mcp_bitscore else "NA",
                    pfam_id=o.best_pfam_acc or "NA",
                    pfam_name=o.best_pfam_name or "NA",
                    pfam_bitscore=f"{o.best_pfam_bitscore:.1f}" if o.best_pfam_bitscore else "NA",
                    pfam_evalue=f"{o.best_pfam_evalue:.2e}" if math.isfinite(o.best_pfam_evalue) else "NA",
                ))


def write_markerout(plvs: List[Plv], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "plv_name", "feature", "name", "start", "end", "strand", "e_value", "score"])
        for p in plvs:
            writer.writerow([p.contig, p.plv_id, "PLV", ".", p.start, p.end, ".", "NA", p.length])
            if p.tir:
                writer.writerow([p.contig, p.plv_id, "TIR_left", ".", p.tir.left_start, p.tir.left_end, "+", p.tir.tir_evalue, p.tir.score])
                writer.writerow([p.contig, p.plv_id, "TIR_right", ".", p.tir.right_start, p.tir.right_end, "-", p.tir.tir_evalue, p.tir.score])
            if p.tsd and p.tir:
                l_end = p.tir.left_start - 1 - p.tsd.left_shift
                l_start = l_end - p.tsd.length + 1
                r_start = p.tir.right_end + 1 + p.tsd.right_shift
                r_end = r_start + p.tsd.length - 1
                writer.writerow([p.contig, p.plv_id, "TSD_5p", p.tsd.sequence_left, l_start, l_end, "+", "NA", f"{p.tsd.identity:.1f}"])
                writer.writerow([p.contig, p.plv_id, "TSD_3p", p.tsd.sequence_right, r_start, r_end, "+", "NA", f"{p.tsd.identity:.1f}"])
            for o in p.orfs:
                feature = "marker" if _is_plv_hallmark(o) else "orf"
                name = _annotation_name(o)
                ev = o.family_evalue if _is_plv_hallmark(o) else o.best_pfam_evalue
                score = o.family_bitscore if _is_plv_hallmark(o) else o.best_pfam_bitscore
                writer.writerow([
                    p.contig, p.plv_id, feature, name, o.start, o.end,
                    "+" if o.strand >= 0 else "-",
                    f"{ev:.2e}" if math.isfinite(ev) else "NA",
                    f"{score:.1f}" if score else "NA",
                ])


def write_nucleotide_fasta(plvs: List[Plv], fa: pyfastx.Fasta, genome_path: Path, path: Path) -> None:
    with open(path, "w") as fh:
        for p in plvs:
            seq = _fetch_seq(fa, genome_path, p.contig, p.start, p.end)
            fh.write(
                f">{p.plv_id} contig={p.contig} start={p.start} end={p.end} "
                f"length={p.length} boundary={p.boundary_method} has_tir={'yes' if p.has_tir else 'no'}\n"
            )
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def write_protein_fasta(plvs: List[Plv], path: Path) -> None:
    with open(path, "w") as fh:
        for p in plvs:
            for i, o in enumerate(p.orfs, start=1):
                label = f"orf{i:05d}"
                annotation = _annotation_name(o)
                marker = f" marker={annotation}" if annotation else ""
                fh.write(f">{p.plv_id}_{label}{marker} length={len(o.protein)}\n")
                for j in range(0, len(o.protein), 80):
                    fh.write(o.protein[j:j + 80] + "\n")


def write_cds_fasta(plvs: List[Plv], fa: pyfastx.Fasta, genome_path: Path, path: Path) -> None:
    with open(path, "w") as fh:
        for p in plvs:
            for i, o in enumerate(p.orfs, start=1):
                label = f"orf{i:05d}"
                seq = _fetch_seq(fa, genome_path, o.contig, o.start, o.end)
                if o.strand < 0:
                    seq = _revcomp(seq)
                annotation = _annotation_name(o)
                marker = f" marker={annotation}" if annotation else ""
                fh.write(f">{p.plv_id}_{label}{marker} length={len(seq)}\n")
                for j in range(0, len(seq), 80):
                    fh.write(seq[j:j + 80] + "\n")


def write_marker_peps(plvs: List[Plv], outdir: Path, prefix: str) -> List[Path]:
    marker_dir = outdir / "marker"
    marker_dir.mkdir(parents=True, exist_ok=True)
    by_marker: Dict[str, Dict[str, Orf]] = defaultdict(dict)
    for p in plvs:
        for o in p.orfs:
            if not _is_plv_hallmark(o):
                continue
            current = by_marker[o.family].get(p.plv_id)
            if current is None or o.family_bitscore > current.family_bitscore:
                by_marker[o.family][p.plv_id] = o
    written = []
    for marker in sorted(by_marker, key=_hallmark_sort_key):
        path = marker_dir / f"{prefix}.{marker}.pep"
        with open(path, "w") as fh:
            for plv_id, o in sorted(by_marker[marker].items(), key=lambda kv: _natural_key(kv[0])):
                fh.write(f">{plv_id}_{marker}\n")
                for i in range(0, len(o.protein), 80):
                    fh.write(o.protein[i:i + 80] + "\n")
        written.append(path)
    return written


def _gff_escape(value: str) -> str:
    """Escape a GFF3 attribute value without introducing separators."""
    from urllib.parse import quote
    return quote(str(value), safe="._:-|+*")


def write_gff3(plvs: List[Plv], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for p in plvs:
            pid = _gff_escape(p.plv_id)
            fh.write(
                f"{p.contig}\tfindPLV\tPLV\t{p.start}\t{p.end}\t.\t+\t.\t"
                f"ID={pid};Name=PLV;boundary={_gff_escape(p.boundary_method)};"
                f"has_tir={'yes' if p.has_tir else 'no'}\n"
            )

            if p.tir:
                fh.write(
                    f"{p.contig}\tfindPLV\tTIR\t{p.tir.left_start}\t{p.tir.left_end}\t"
                    f"{p.tir.score:.1f}\t+\t.\tID={pid}.TIR_left;Parent={pid};Name=TIR\n"
                )
                fh.write(
                    f"{p.contig}\tfindPLV\tTIR\t{p.tir.right_start}\t{p.tir.right_end}\t"
                    f"{p.tir.score:.1f}\t-\t.\tID={pid}.TIR_right;Parent={pid};Name=TIR\n"
                )

            if p.tsd and p.tir:
                l_end = p.tir.left_start - 1 - p.tsd.left_shift
                l_start = l_end - p.tsd.length + 1
                r_start = p.tir.right_end + 1 + p.tsd.right_shift
                r_end = r_start + p.tsd.length - 1
                fh.write(
                    f"{p.contig}\tfindPLV\tTSD\t{l_start}\t{l_end}\t.\t+\t.\t"
                    f"ID={pid}.TSD_5p;Parent={pid};Name=TSD;"
                    f"sequence={_gff_escape(p.tsd.sequence_left)}\n"
                )
                fh.write(
                    f"{p.contig}\tfindPLV\tTSD\t{r_start}\t{r_end}\t.\t+\t.\t"
                    f"ID={pid}.TSD_3p;Parent={pid};Name=TSD;"
                    f"sequence={_gff_escape(p.tsd.sequence_right)}\n"
                )

            for i, o in enumerate(p.orfs, start=1):
                orf_label = f"orf{i:05d}"
                name = _annotation_name(o)
                if _is_plv_hallmark(o):
                    score = o.family_bitscore
                    score_field = f"{score:.1f}"
                    attrs = (
                        f"ID={pid}.{orf_label};Parent={pid};"
                        f"Name={_gff_escape(name)};Orf={orf_label};"
                        f"hallmark={_gff_escape(o.family)}"
                    )
                elif o.best_pfam_name:
                    score_field = f"{o.best_pfam_bitscore:.1f}"
                    attrs = (
                        f"ID={pid}.{orf_label};Parent={pid};"
                        f"Name={_gff_escape(name)};Orf={orf_label};"
                        f"Pfam={_gff_escape(o.best_pfam_acc or o.best_pfam_name)}"
                    )
                else:
                    score_field = "."
                    attrs = (
                        f"ID={pid}.{orf_label};Parent={pid};"
                        f"Name=hypothetical;Orf={orf_label}"
                    )
                fh.write(
                    f"{p.contig}\tfindPLV\tCDS\t{o.start}\t{o.end}\t"
                    f"{score_field}\t{'+' if o.strand >= 0 else '-'}\t0\t{attrs}\n"
                )


def log_result_summary(plvs: List[Plv], genome_path: Path) -> None:
    _LOG.info("Result Summary")
    _LOG.info(f"Input genome: {genome_path}")
    _LOG.info(f"PLV candidates: {len(plvs)}")
    if not plvs:
        return
    n_tir = sum(p.has_tir for p in plvs)
    _LOG.info(f"  With TIR:    {n_tir} ({100 * n_tir / len(plvs):.1f}%)")
    _LOG.info(f"  Without TIR: {len(plvs) - n_tir} ({100 * (len(plvs) - n_tir) / len(plvs):.1f}%)")
    lengths = [p.length for p in plvs]
    _LOG.info(f"PLV length: min={min(lengths):,}; max={max(lengths):,}; mean={sum(lengths) / len(lengths):,.0f} bp")
    confidence_counts = Counter(calculate_confidence(p)[1] for p in plvs)
    _LOG.info("Confidence: " + ", ".join(f"{k}={confidence_counts[k]}" for k in ("HIGH", "MEDIUM", "LOW") if confidence_counts[k]))
    fams = Counter(f for p in plvs for f in p.families_present)
    if fams:
        _LOG.info("Most frequent markers: " + ", ".join(f"{k}({v})" for k, v in fams.most_common(6)))


def _write_no_plv_notice(outdir: Path, reason: str) -> None:
    path = outdir / "No_PLV_was_found.txt"
    path.write_text(
        "No PLV was found in this run.\n"
        "Please read run.log in this output directory for full details.\n"
        f"Reason: {reason}\n"
    )
    _LOG.output(f"No-PLV notice -> {path}")


def _get_tool_versions() -> Dict[str, str]:
    versions: Dict[str, str] = {}
    for pkg in ("pyfastx", "pyhmmer", "pyrodigal"):
        try:
            versions[pkg] = _pkg_version(pkg)
        except PackageNotFoundError:
            versions[pkg] = "unknown"
    try:
        proc = subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=False)
        lines = (proc.stdout or proc.stderr).strip().splitlines()
        versions["blastn"] = lines[0] if lines else "unknown"
    except OSError:
        versions["blastn"] = "not_found"
    versions["python"] = sys.version.split()[0]
    return versions

class _Parser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HELP_TEXT

    def format_usage(self) -> str:
        return USAGE_TEXT


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = _Parser(prog="findPLV_v3.py", add_help=False)
    p.add_argument("-h", "--help", action="help")
    p.add_argument("genome", type=Path)
    p.add_argument("-db", "--db", type=Path, required=True)
    p.add_argument("--prefix", type=str, required=True)
    p.add_argument("-o", "--outdir", type=Path, default=Path(f"Result_{datetime.now().strftime('%Y%m%d')}"))
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("-p", "--parallel", type=int, default=None)
    p.add_argument("-g", "--gff", type=Path, default=None)
    p.add_argument("-e", "--evalue", type=float, default=DEFAULTS["hmm_evalue"])
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    cfg = dict(DEFAULTS)
    cfg["hmm_evalue"] = args.evalue

    if args.threads < 1:
        print("--threads must be >= 1", file=sys.stderr)
        return 2
    parallel = max(1, int(args.parallel if args.parallel is not None else args.threads))

    args.outdir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.outdir / "run.log")
    t0 = time.time()

    _LOG.info(
        f"findPLV_v3 started | prefix='{args.prefix}' | threads={args.threads} | "
        f"parallel={parallel} | genome={args.genome}"
    )
    _LOG.info(
        f"Parameters | min_contig={cfg['min_contig']:,} | min_plv_length={cfg['min_plv_length']:,} | "
        f"max_plv_length={cfg['max_plv_length']:,} | seed_window={cfg['seed_window']:,} | "
        f"min_families={cfg['min_families']} | TIR={cfg['tir_min_len']}-{cfg['tir_max_len']} bp, "
        f">={cfg['tir_min_id']:.1f}% id | HMM E={cfg['hmm_evalue']:g}"
    )

    if not args.genome.is_file():
        _LOG.error(f"Genome file not found: {args.genome}")
        return 2
    if not args.db.is_dir():
        _LOG.error(f"Database directory not found: {args.db}")
        return 2
    if args.gff is not None and not args.gff.is_file():
        _LOG.error(f"GFF file not found: {args.gff}")
        return 2

    db_paths = {key: args.db / filename for key, filename in DB_FILES.items()}
    missing = [path for path in db_paths.values() if not path.is_file()]
    if missing:
        for path in missing:
            _LOG.error(f"Required database file not found: {path}")
        return 2
    if not shutil.which("blastn"):
        _LOG.error("blastn not found in PATH. Install BLAST+ (e.g. conda install -c bioconda blast)")
        return 2

    versions = _get_tool_versions()
    _LOG.info("Tool versions | " + " | ".join(f"{k}={v}" for k, v in versions.items()))
    _LOG.info(f"Command line | {' '.join(sys.argv)}")

    try:
        orfs_by_id, contig_lengths = predict_orfs(args.genome, cfg["min_contig"], args.threads)
    except Exception as exc:
        _LOG.error(str(exc))
        return 2

    host_intervals: Dict[str, List[Tuple[int, int]]] = {}
    if args.gff is not None:
        host_intervals = parse_host_gff_intervals(args.gff)
        filter_orfs_by_host_gff(orfs_by_id, host_intervals)
        if not orfs_by_id:
            _write_no_plv_notice(args.outdir, "Host GFF mask removed all predicted ORFs")
            log_result_summary([], args.genome)
            return 0

    try:
        scan_hallmarks(orfs_by_id, db_paths["hallmarks"], cfg["hmm_evalue"], args.threads)
        scan_nonplv_mcp(orfs_by_id, db_paths["mcp_nonplv"], cfg["hmm_evalue"], args.threads)
        n_mcp, _n_amb, _n_other = classify_hallmarks(orfs_by_id, cfg["mcp_score_margin"])
    except Exception as exc:
        _LOG.error(f"PLV hallmark/MCP HMM scanning failed: {exc}")
        return 2

    if n_mcp == 0:
        _write_no_plv_notice(args.outdir, "No competitively classified PLV MCP was detected")
        log_result_summary([], args.genome)
        return 0

    clusters = find_seed_regions(
        orfs_by_id, cfg["seed_window"], cfg["min_families"],
        cfg["cluster_merge_gap"], cfg["max_cluster_span"],
    )
    if not clusters:
        _write_no_plv_notice(args.outdir, "No PLV marker cluster passed the seeding criteria")
        log_result_summary([], args.genome)
        return 0

    orfs_by_contig: Dict[str, List[Orf]] = defaultdict(list)
    for o in orfs_by_id.values():
        orfs_by_contig[o.contig].append(o)
    for orfs in orfs_by_contig.values():
        orfs.sort(key=lambda o: o.start)

    pfam_target_ids = set()
    pfam_flank = cfg["max_plv_length"]
    for cluster in clusters:
        left = max(1, cluster["cluster_start"] - pfam_flank)
        right = cluster["cluster_end"] + pfam_flank
        for o in orfs_by_contig.get(cluster["contig"], []):
            if o.start > right:
                break
            if o.end >= left:
                pfam_target_ids.add(o.orf_id)
    pfam_targets = [o for o in orfs_by_id.values() if o.orf_id in pfam_target_ids]
    try:
        scan_pfam(pfam_targets, db_paths["pfam"], cfg["hmm_evalue"], args.threads)
    except Exception as exc:
        _LOG.error(f"Pfam-A scanning failed: {exc}")
        return 2

    tasks = [
        dict(
            cluster_index=i, contig=cl["contig"], cluster_start=cl["cluster_start"],
            cluster_end=cl["cluster_end"], contig_length=contig_lengths[cl["contig"]],
            cfg=cfg, genome_path=str(args.genome), contig_orfs=orfs_by_contig[cl["contig"]],
            host_intervals=host_intervals.get(cl["contig"], []),
        )
        for i, cl in enumerate(clusters, start=1)
    ]

    n_workers = max(1, min(parallel, len(tasks)))
    _LOG.info(f"Candidate processing: {len(tasks)} cluster(s) across {n_workers} worker(s); BLASTN uses 1 thread per worker")
    executor = None
    if n_workers > 1:
        executor = ProcessPoolExecutor(max_workers=n_workers)
        results_iter = executor.map(_process_cluster, tasks, chunksize=1)
    else:
        results_iter = (_process_cluster(t) for t in tasks)

    plvs: List[Plv] = []
    try:
        for res in results_iter:
            for level, msg in res.get("log_msgs", []):
                (_LOG.warning if level == "warning" else _LOG.info)(msg)
            if res.get("status") == "ok":
                p = res["plv"]
                plvs.append(p)
                _LOG.info(
                    f"Cluster {res['cluster_index']} accepted | {p.contig}:{p.start:,}-{p.end:,} | "
                    f"{p.length:,} bp | markers={p.n_families} | boundary={p.boundary_method} | "
                    f"TIR={'yes' if p.has_tir else 'no'} | TSD={_tsd_status(p.tsd)}"
                )
            else:
                _LOG.warning(res.get("message", f"Cluster {res.get('cluster_index', '?')} discarded"))
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    plvs = deduplicate_plvs(plvs, cfg["dedup_min_reciprocal_overlap"])
    plvs.sort(key=lambda p: (p.contig, p.start, p.end))
    for i, p in enumerate(plvs, start=1):
        p.plv_id = f"{args.prefix}_PLV_{i:03d}"

    if not plvs:
        _write_no_plv_notice(args.outdir, "All candidate clusters failed final PLV filters")
        log_result_summary([], args.genome)
        return 0

    fa = pyfastx.Fasta(str(args.genome), build_index=True, uppercase=True)
    summary_path = args.outdir / f"{args.prefix}.summary.tsv"
    func_path = args.outdir / f"{args.prefix}.func.tsv"
    markerout_path = args.outdir / f"{args.prefix}.markerout"
    fna_path = args.outdir / f"{args.prefix}.plv.fna"
    pep_path = args.outdir / f"{args.prefix}.plv.pep"
    cds_path = args.outdir / f"{args.prefix}.plv.cds"
    gff_path = args.outdir / f"{args.prefix}.plv.gff3"

    write_summary_tsv(plvs, summary_path)
    write_func_tsv(plvs, func_path)
    write_markerout(plvs, markerout_path)
    write_nucleotide_fasta(plvs, fa, args.genome, fna_path)
    write_protein_fasta(plvs, pep_path)
    write_cds_fasta(plvs, fa, args.genome, cds_path)
    marker_paths = write_marker_peps(plvs, args.outdir, args.prefix)
    write_gff3(plvs, gff_path)

    log_result_summary(plvs, args.genome)
    _LOG.output(f"Summary table -> {summary_path}")
    _LOG.output(f"Function table -> {func_path}")
    _LOG.output(f"Marker table  -> {markerout_path}")
    _LOG.output(f"PLV FASTA     -> {fna_path}")
    _LOG.output(f"Protein FASTA -> {pep_path}")
    _LOG.output(f"CDS FASTA     -> {cds_path}")
    _LOG.output(f"GFF3          -> {gff_path}")
    _LOG.output(f"Marker FASTAs -> {args.outdir / 'marker'} ({len(marker_paths)} file(s))")
    _LOG.output(f"Run log       -> {args.outdir / 'run.log'}")
    _LOG.info(f"findPLV_v3 completed in {time.time() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
