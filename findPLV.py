#!/usr/bin/env python3
"""findPLV.py — Identify Polinton-Like Viruses in eukaryotic genome assemblies.

Usage: findPLV.py genome.fa --mcp-plv MCP_PLV.hmm --mcp-ncldv-virophage NCLDV.hmm
                             --pfam Pfam-A.hmm [OPTIONS]
Author:
  Dede Kurniawan (dedekurniawan@genomics.cn)
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import numpy as np
import pandas as pd
import pyfastx
import pyhmmer
import pyrodigal


# ── Runtime defaults ──────────────────────────────────────────────────────────

DEFAULTS = dict(
    min_contig                  = 5_000,
    min_plv_length              = 5_000,
    # Seeding
    window_size                 = 20_000,   
    window_step                 = 5_000,    
    min_families                = 3,        
    require_mcp                 = True,
    cluster_merge_gap           = 8_000,    
    max_cluster_span            = 60_000,   
    # TIR detection
    flank_for_tir               = 50_000,
    tir_min_insert              = 6_000,
    tir_max_insert              = 80_000,
    tir_min_len                 = 10,
    tir_max_len                 = 8_000,
    tir_min_id                  = 0.0,
    # TSD
    tsd_min                     = 4,
    tsd_max                     = 12,
    tsd_max_slide               = 2,
    # GC
    gc_window                   = 500,
    gc_flank                    = 50_000,
    # E-value thresholds
    mcp_evalue                  = 1e-5,
    mcp_minor_evalue            = 1e-3,
    pfam_evalue                 = 1e-5,
    ncldv_virophage_mcp_evalue  = 1e-5,
    # einverted parameters
    einverted_match             = 3,
    einverted_mismatch          = -4,
    einverted_gap               = 12,
    einverted_threshold         = 15,
    einverted_maxrepeat         = 80_000,
    # TIR-less span control (anchor-driven, three-phase, GC-aware)
    max_plv_length              = 35_000,   
    edge_gap_trim               = 5_000,    
    gc_flank_search             = 10_000,   
                                            
    gc_min_delta_pct            = 3.0,      
                                            
    # Deduplication
    dedup_min_reciprocal_overlap = 0.50,    
)

# ── PLV marker allow-list (Pfam accession → label, counts_toward_seeding) ────

DEFAULT_PLV_FAMILIES: List[Tuple[str, str, bool]] = [
    ("PF03175",   "pPolB",            True),
    ("PF00665",   "RVE-INT",          True),
    ("PF04665",   "A32-pATPase",      True),
    ("PF00589",   "Y-rec",            True),
    ("PF19835",   "GIY-YIG",          True),
    ("PF04735",   "TVpol-S3H",        True),
    ("PF13385",   "Laminin_G_3",      True),
    ("PF03903",   "Peptidase_S74",    True),
    ("PF02902",   "vUlp1-PRO",        True),
    ("PF13461",   "AEP-primase",      True),
    ("PF00271",   "S1H-helicase",     True),
    ("PF00476",   "PolA",             True),
]

MCP_LABEL       = "MCP"   
MCP_MINOR_LABEL = "mCP"   

ANCHOR_FAMILIES_BASE: frozenset = frozenset({
    MCP_LABEL,        
    MCP_MINOR_LABEL, 
    "pPolB",         
    "A32-pATPase",   
    "RVE-INT",       
})


SUPPORTIVE_FAMILIES: frozenset = frozenset({
    "Y-rec",          # tyrosine recombinase
    "GIY-YIG",        # GIY-YIG endonuclease
    "TVpol-S3H",      # TVpol superfamily 3 helicase
    "Laminin_G_3",    # laminin G domain
    "Peptidase_S74",  # S74 protease
    "vUlp1-PRO",      # viral Ulp1-like protease
    "AEP-primase",    # archaeo-eukaryotic primase
    "S1H-helicase",   # SF1H helicase
    "PolA",           # family-A DNA polymerase
})

_CAPSID_LABELS: frozenset = frozenset({MCP_LABEL, MCP_MINOR_LABEL})


# ── Logging helpers ───────────────────────────────────────────────────────────
_LOG = logging.getLogger("findPLV")


def setup_logging(log_path: Optional[Path] = None) -> None:
    _LOG.setLevel(logging.INFO)
    fmt = logging.Formatter("%(message)s")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    _LOG.addHandler(sh)
    if log_path is not None:
        fh = logging.FileHandler(log_path, mode="w")
        fh.setFormatter(fmt)
        _LOG.addHandler(fh)


def log_info(msg: str) -> None:    _LOG.info(f"[Info] {msg}")
def log_output(msg: str) -> None:  _LOG.info(f"[Output] {msg}")
def log_warning(msg: str) -> None: _LOG.warning(f"[Warning] {msg}")
def log_error(msg: str) -> None:   _LOG.error(f"[Error] {msg}")


# ── Data-classes ──────────────────────────────────────────────────────────────
@dataclass
class Orf:
    orf_id: str
    contig: str
    start: int           # 1-based inclusive
    end: int             # 1-based inclusive
    strand: int          # +1 / -1
    protein: str
    partial: bool
    family: Optional[str] = None
    family_bitscore: float = 0.0
    family_evalue: float = float("inf")
    family_source: Optional[str] = None   # "MCP" | "mCP" | "Pfam"
    pfam_hits: List[Tuple[str, str, float, float]] = field(default_factory=list)
    best_pfam_acc: Optional[str] = None
    best_pfam_label: Optional[str] = None
    best_pfam_bitscore: float = 0.0
    best_pfam_evalue: float = float("inf")


@dataclass
class TirPair:
    left_start: int   # genome-relative, 1-based inclusive
    left_end: int
    right_start: int
    right_end: int
    tir_length: int
    insert_size: int
    plv_span: int
    tir_identity: float
    score: int
    matches: int
    total: int
    gaps: int
    left_seq: str = ""
    right_seq: str = ""


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
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    length: int
    tir: Optional[TirPair]
    tsd: Optional[Tsd]
    orfs: List[Orf]
    n_families: int
    families_present: List[str]
    mcp_best_bitscore: float
    gc_plv: float
    has_tir: bool = False


# ── Stage 1: ORF prediction ───────────────────────────────────────────────────
def _predict_orfs_on_contig(
    args: Tuple[str, str],
) -> Tuple[str, List[Tuple[str, int, int, int, str, bool]], Optional[str]]:
    contig_name, seq = args
    try:
        gene_finder = pyrodigal.GeneFinder(meta=True)
        genes = gene_finder.find_genes(seq.encode("ascii"))
    except Exception as e:
        return contig_name, [], f"pyrodigal failed on {contig_name}: {e}"
    out: List[Tuple[str, int, int, int, str, bool]] = []
    for i, gene in enumerate(genes, start=1):
        protein = gene.translate().rstrip("*")
        if not protein:
            continue
        partial = bool(
            getattr(gene, "partial_begin", False) or getattr(gene, "partial_end", False)
        )
        out.append((f"orf{i:05d}", int(gene.begin), int(gene.end),
                    int(gene.strand), protein, partial))
    return contig_name, out, None


def predict_orfs(
    genome_path: Path, min_contig: int, threads: int
) -> Tuple[Dict[str, Orf], Dict[str, int]]:
    fa = pyfastx.Fasta(str(genome_path), build_index=True, uppercase=True)
    orfs_by_id: Dict[str, Orf] = {}
    contig_lengths: Dict[str, int] = {}
    n_kept = n_skipped = 0
    work_items: List[Tuple[str, str]] = []

    for record in fa:
        cname = record.name
        seqlen = len(record.seq)
        contig_lengths[cname] = seqlen
        if seqlen < min_contig:
            n_skipped += 1
            continue
        n_kept += 1
        work_items.append((cname, str(record.seq)))

    if not work_items:
        log_error("No contigs passed the minimum-length filter. Check input FASTA and --min-contig.")
        sys.exit(1)

    if threads <= 1 or len(work_items) == 1:
        results_iter = (_predict_orfs_on_contig(wi) for wi in work_items)
    else:
        executor = ProcessPoolExecutor(max_workers=threads)
        results_iter = executor.map(_predict_orfs_on_contig, work_items, chunksize=1)

    total_orfs = 0
    try:
        for cname, orf_records, err in results_iter:
            if err:
                log_warning(f"{err}; skipping contig")
                continue
            for suffix, start, end, strand, protein, partial in orf_records:
                orf_id = f"{cname}__{suffix}"
                orfs_by_id[orf_id] = Orf(
                    orf_id=orf_id, contig=cname,
                    start=start, end=end, strand=strand,
                    protein=protein, partial=partial,
                )
                total_orfs += 1
    finally:
        if threads > 1 and len(work_items) > 1:
            executor.shutdown(wait=True)

    log_info(
        f"ORF prediction: {total_orfs:,} ORFs on {n_kept:,} contig(s) "
        f"(>= {min_contig:,} bp; {n_skipped:,} skipped)"
    )
    if total_orfs == 0:
        log_error("No ORFs predicted. Check input FASTA and --min-contig.")
        sys.exit(1)
    return orfs_by_id, contig_lengths


# ── Stage 2: HMM scanning ─────────────────────────────────────────────────────
def _digital_sequences(
    orfs: Iterable[Orf], alphabet: pyhmmer.easel.Alphabet
) -> List[pyhmmer.easel.DigitalSequence]:
    seqs = []
    for o in orfs:
        try:
            t = pyhmmer.easel.TextSequence(name=o.orf_id.encode(), sequence=o.protein)
            seqs.append(t.digitize(alphabet))
        except Exception as e:
            log_warning(f"Skipping ORF {o.orf_id} during digitization: {e}")
    return seqs


def _query_meta(top_hits) -> Tuple[str, Optional[str]]:
    try:
        name_attr = top_hits.query.name
        name = name_attr.decode() if isinstance(name_attr, bytes) else name_attr
        acc_attr = top_hits.query.accession
        acc = (acc_attr.decode() if isinstance(acc_attr, bytes) else acc_attr) if acc_attr else None
    except AttributeError:
        name_attr = top_hits.query_name
        name = name_attr.decode() if isinstance(name_attr, bytes) else name_attr
        acc_attr = top_hits.query_accession
        acc = (acc_attr.decode() if isinstance(acc_attr, bytes) else acc_attr) if acc_attr else None
    if acc:
        acc = acc.split(".")[0]
    return name, acc


def scan_mcp(
    orfs_by_id: Dict[str, Orf], mcp_hmm_path: Path, evalue: float, threads: int
) -> None:
    """Scan all proteins against the major capsid (MCP) HMM; annotates Orf.family in place."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_sequences(orfs_by_id.values(), alphabet)
    with pyhmmer.plan7.HMMFile(str(mcp_hmm_path)) as hf:
        hmms = list(hf)
    if not hmms:
        log_error(f"No HMM profile found in {mcp_hmm_path}")
        sys.exit(1)
    n_hits = 0
    for top_hits in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
        for hit in top_hits:
            if not hit.included:
                continue
            target = hit.name
            target = target.decode() if isinstance(target, bytes) else target
            o = orfs_by_id.get(target)
            if o is None:
                continue
            if hit.score > o.family_bitscore:
                o.family         = MCP_LABEL
                o.family_bitscore = float(hit.score)
                o.family_evalue  = float(hit.evalue)
                o.family_source  = "MCP"
            n_hits += 1
    log_info(f"MCP HMM scan: {n_hits:,} hit(s) at E-value <= {evalue:g}")


def scan_mcp_minor(
    orfs_by_id: Dict[str, Orf],
    mcp_minor_hmm_path: Path,
    evalue: float,
    threads: int,
) -> None:
    """Scan all proteins against the minor capsid (mCP) HMM; annotates Orf.family in place.

    An ORF already assigned MCP retains that label regardless of bitscore
    (MCP always takes priority over mCP). For all other ORFs, mCP wins on
    the usual bitscore comparison. Run this AFTER scan_mcp().
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_sequences(orfs_by_id.values(), alphabet)
    with pyhmmer.plan7.HMMFile(str(mcp_minor_hmm_path)) as hf:
        hmms = list(hf)
    if not hmms:
        log_error(f"No HMM profile found in {mcp_minor_hmm_path}")
        sys.exit(1)
    n_hits = 0
    for top_hits in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
        for hit in top_hits:
            if not hit.included:
                continue
            target = hit.name
            target = target.decode() if isinstance(target, bytes) else target
            o = orfs_by_id.get(target)
            if o is None:
                continue
            # MCP takes priority: never let mCP overwrite a confirmed major capsid.
            if o.family == MCP_LABEL:
                continue
            if hit.score > o.family_bitscore:
                o.family          = MCP_MINOR_LABEL
                o.family_bitscore = float(hit.score)
                o.family_evalue   = float(hit.evalue)
                o.family_source   = "mCP"
                n_hits += 1
    log_info(f"mCP HMM scan: {n_hits:,} hit(s) at E-value <= {evalue:g}")


def scan_pfam(
    orfs_to_scan: List[Orf],
    pfam_hmm_path: Path,
    family_table: Dict[str, Tuple[str, bool]],
    evalue: float,
    threads: int,
) -> None:
    """Scan a subset of proteins against Pfam-A.

    Records all hits to Orf.pfam_hits; assigns Orf.family only when the
    bitscore exceeds the current best AND the ORF has not already been
    assigned a capsid label (MCP or mCP).
    """
    if not orfs_to_scan:
        return
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_sequences(orfs_to_scan, alphabet)
    by_id = {o.orf_id: o for o in orfs_to_scan}

    n_hits_total = n_hits_plv = 0
    with pyhmmer.plan7.HMMFile(str(pfam_hmm_path)) as hf:
        for top_hits in pyhmmer.hmmsearch(hf, seqs, cpus=threads, E=evalue):
            hmm_name, hmm_acc = _query_meta(top_hits)
            for hit in top_hits:
                if not hit.included:
                    continue
                target = hit.name
                target = target.decode() if isinstance(target, bytes) else target
                o = by_id.get(target)
                if o is None:
                    continue

                acc_rec = hmm_acc if hmm_acc else hmm_name
                if hmm_acc and hmm_acc in family_table:
                    label_rec = family_table[hmm_acc][0]
                else:
                    label_rec = hmm_name

                o.pfam_hits.append(
                    (acc_rec, label_rec, float(hit.score), float(hit.evalue))
                )
                n_hits_total += 1

                if float(hit.score) > o.best_pfam_bitscore:
                    o.best_pfam_bitscore = float(hit.score)
                    o.best_pfam_evalue   = float(hit.evalue)
                    o.best_pfam_acc      = acc_rec
                    o.best_pfam_label    = label_rec

                if hmm_acc and hmm_acc in family_table:
                    n_hits_plv += 1

                # Capsid labels (MCP / mCP) are protected from Pfam overwriting.
                if hit.score > o.family_bitscore and o.family not in _CAPSID_LABELS:
                    o.family         = label_rec
                    o.family_bitscore = float(hit.score)
                    o.family_evalue  = float(hit.evalue)
                    o.family_source  = "Pfam"

    log_info(
        f"Pfam-A scan: {n_hits_total:,} total hit(s) "
        f"({n_hits_plv:,} on curated PLV markers; "
        f"{n_hits_total - n_hits_plv:,} additional hits retained)"
    )


# ── Stage 2c: multi-anchor seeding with two-pass merge ────────────────────────
def find_seed_regions(
    orfs_by_id: Dict[str, Orf],
    anchor_families: frozenset,
    window_size: int,
    min_families: int,
    require_mcp: bool,
    cluster_merge_gap: int,
    max_cluster_span: int,
) -> List[dict]:
    """Identify candidate PLV marker clusters via sliding windows.

    Algorithm
    ---------
    For every anchor-family ORF (MCP, mCP, pPolB, A32-pATPase, RVE-INT)
    on each contig, centre a window of `window_size` bp.  Collect all
    labeled ORFs inside the window.  Keep the window when:
      - it contains ≥ min_families distinct families, AND
      - (if require_mcp) at least one MCP ORF is present.

    Pass 1 — merge overlapping raw windows on the same contig.

    Pass 2 — merge adjacent clusters within `cluster_merge_gap` bp of
    each other, provided the merged span ≤ `max_cluster_span`.

    After both passes, any cluster lacking MCP is dropped (if require_mcp).
    """
    by_contig: Dict[str, List[Orf]] = {}
    for o in orfs_by_id.values():
        by_contig.setdefault(o.contig, []).append(o)
    for v in by_contig.values():
        v.sort(key=lambda x: x.start)

    half = window_size // 2
    raw_windows: List[dict] = []

    for contig, orfs in by_contig.items():
        labeled = [o for o in orfs if o.family is not None]
        if not labeled:
            continue
        # Only anchor-family ORFs generate windows.
        anchors = [o for o in labeled if o.family in anchor_families]
        if not anchors:
            continue

        l_starts = np.array([o.start for o in labeled])
        l_ends   = np.array([o.end   for o in labeled])

        for a in anchors:
            mid    = (a.start + a.end) // 2
            wstart = mid - half
            wend   = mid + half
            mask   = (l_starts <= wend) & (l_ends >= wstart)
            members = [labeled[i] for i in np.nonzero(mask)[0]]
            fams    = {m.family for m in members}

            if len(fams) < min_families:
                continue
            if require_mcp and MCP_LABEL not in fams:
                continue

            raw_windows.append(dict(
                contig        = contig,
                cluster_start = min(m.start for m in members),
                cluster_end   = max(m.end   for m in members),
                orf_ids       = [m.orf_id for m in members],
                families      = sorted(fams),
            ))

    # ── Pass 1: collapse overlapping windows ──────────────────────────────────
    raw_windows.sort(key=lambda c: (c["contig"], c["cluster_start"]))
    pass1: List[dict] = []
    for w in raw_windows:
        if (pass1
                and pass1[-1]["contig"] == w["contig"]
                and w["cluster_start"] <= pass1[-1]["cluster_end"]):
            m = pass1[-1]
            m["cluster_end"] = max(m["cluster_end"], w["cluster_end"])
            m["orf_ids"]     = sorted(set(m["orf_ids"]) | set(w["orf_ids"]))
            m["families"]    = sorted(set(m["families"]) | set(w["families"]))
        else:
            pass1.append(dict(w))

    # ── Pass 2: merge clusters within cluster_merge_gap ───────────────────────
    merged: List[dict] = []
    for c in pass1:
        if merged and merged[-1]["contig"] == c["contig"]:
            gap = c["cluster_start"] - merged[-1]["cluster_end"] - 1
            prospective_span = c["cluster_end"] - merged[-1]["cluster_start"] + 1
            if gap <= cluster_merge_gap and prospective_span <= max_cluster_span:
                m = merged[-1]
                m["cluster_end"] = max(m["cluster_end"], c["cluster_end"])
                m["orf_ids"]     = sorted(set(m["orf_ids"]) | set(c["orf_ids"]))
                m["families"]    = sorted(set(m["families"]) | set(c["families"]))
                continue
        merged.append(dict(c))

    # ── Final MCP check after all merging ─────────────────────────────────────
    if require_mcp:
        before = len(merged)
        merged = [c for c in merged if MCP_LABEL in c["families"]]
        n_dropped = before - len(merged)
        if n_dropped:
            log_info(
                f"Seeding: dropped {n_dropped} merged cluster(s) lacking MCP "
                f"after gap-merge"
            )

    log_info(
        f"Seeding: {len(merged)} candidate PLV region(s) "
        f"[anchors={','.join(sorted(anchor_families))} | "
        f"window={window_size:,} bp | "
        f"min_families={min_families} | "
        f"merge_gap={cluster_merge_gap:,} bp | "
        f"max_span={max_cluster_span:,} bp]"
    )
    return merged


# ── Stage 2d: NCLDV / Virophage MCP exclusion ────────────────────────────────
def scan_ncldv_virophage_mcp(
    orfs: List[Orf],
    hmm_path: Path,
    evalue: float,
    threads: int,
) -> Dict[str, List[Tuple[str, float, float]]]:
    if not orfs:
        return {}
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = _digital_sequences(orfs, alphabet)
    if not seqs:
        return {}
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        hmms = list(hf)
    if not hmms:
        log_warning(f"No HMM profile found in {hmm_path}; skipping NCLDV/Virophage MCP filter")
        return {}
    hits_by_orf: Dict[str, List[Tuple[str, float, float]]] = {}
    for top_hits in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
        hmm_name, _ = _query_meta(top_hits)
        for hit in top_hits:
            if not hit.included:
                continue
            target = hit.name
            target = target.decode() if isinstance(target, bytes) else target
            hits_by_orf.setdefault(target, []).append(
                (hmm_name, float(hit.score), float(hit.evalue))
            )
    return hits_by_orf


def filter_clusters_by_ncldv_virophage(
    clusters: List[dict],
    orfs_by_id: Dict[str, Orf],
    hmm_path: Path,
    evalue: float,
    threads: int,
) -> List[dict]:
    if not clusters:
        return clusters
    orfs_by_contig: Dict[str, List[Orf]] = {}
    for o in orfs_by_id.values():
        orfs_by_contig.setdefault(o.contig, []).append(o)
    for v in orfs_by_contig.values():
        v.sort(key=lambda x: x.start)

    candidate_orfs: Dict[str, Orf] = {}
    cluster_orf_ids: List[List[str]] = []
    for cl in clusters:
        cs, ce = cl["cluster_start"], cl["cluster_end"]
        ids_in: List[str] = []
        for o in orfs_by_contig.get(cl["contig"], []):
            if o.start <= ce and o.end >= cs:
                candidate_orfs[o.orf_id] = o
                ids_in.append(o.orf_id)
        cluster_orf_ids.append(ids_in)

    if not candidate_orfs:
        log_info("NCLDV/Virophage MCP filter: no ORFs overlap any seed cluster; nothing to scan")
        return clusters

    log_info(
        f"NCLDV/Virophage MCP filter: scanning {len(candidate_orfs):,} ORF(s) "
        f"from {len(clusters)} cluster(s) at E-value <= {evalue:g}"
    )
    hits_by_orf = scan_ncldv_virophage_mcp(list(candidate_orfs.values()), hmm_path, evalue, threads)

    if not hits_by_orf:
        log_info("NCLDV/Virophage MCP filter: no hits; all seed clusters retained")
        return clusters

    kept: List[dict] = []
    n_ncldv = n_virophage = n_both = 0
    for cl, ids_in in zip(clusters, cluster_orf_ids):
        matched: set = set()
        best_hit: Optional[Tuple[str, str, float, float]] = None
        for oid in ids_in:
            for hmm_name, score, ev in hits_by_orf.get(oid, ()):
                matched.add(hmm_name)
                if best_hit is None or score > best_hit[2]:
                    best_hit = (oid, hmm_name, score, ev)
        if not matched:
            kept.append(cl)
            continue
        is_ncldv     = any(n.lower() == "mcp" for n in matched)
        is_virophage = any("virophage" in n.lower() for n in matched)
        if is_ncldv and is_virophage:
            n_both += 1
            tag = "NCLDV+Virophage MCP"
        elif is_ncldv:
            n_ncldv += 1
            tag = "NCLDV MCP"
        elif is_virophage:
            n_virophage += 1
            tag = "Virophage MCP"
        else:
            tag = ",".join(sorted(matched))
        if best_hit is not None:
            log_info(
                f"  Dropping cluster {cl['contig']}:{cl['cluster_start']:,}-{cl['cluster_end']:,} "
                f"({tag}; best: {best_hit[0]} -> {best_hit[1]}, "
                f"bitscore={best_hit[2]:.1f}, E={best_hit[3]:.1g})"
            )

    n_dropped = len(clusters) - len(kept)
    log_info(
        f"NCLDV/Virophage MCP filter: dropped {n_dropped} cluster(s) "
        f"(NCLDV: {n_ncldv}, Virophage: {n_virophage}, both: {n_both}); "
        f"{len(kept)} retained"
    )
    return kept


# ── Stage 3: TIR detection with einverted ─────────────────────────────────────
EINV_HEADER_RE = re.compile(
    r"^(?P<n>\S+):\s+Score\s+(?P<score>\d+):\s+"
    r"(?P<matches>\d+)/(?P<total>\d+)\s+\((?P<pct>\d+)%\)\s+matches,\s+"
    r"(?P<gaps>\d+)\s+gaps\s*$"
)
EINV_POS_RE = re.compile(r"^\s*(?P<s>\d+)\s+(?P<seq>.+?)\s+(?P<e>\d+)\s*$")


def parse_einverted(inv_path: Path) -> List[TirPair]:
    if not inv_path.exists():
        return []
    text = inv_path.read_text()
    if not text.strip():
        return []
    pairs: List[TirPair] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        h = EINV_HEADER_RE.match(lines[i])
        if not h:
            i += 1
            continue
        pos_matches: List[re.Match] = []
        j, scanned = i + 1, 0
        while j < len(lines) and len(pos_matches) < 2 and scanned < 6:
            pm = EINV_POS_RE.match(lines[j])
            if pm:
                pos_matches.append(pm)
            j += 1
            scanned += 1
        if len(pos_matches) < 2:
            i = max(j, i + 1)
            continue
        try:
            s1, e1 = int(pos_matches[0].group("s")), int(pos_matches[0].group("e"))
            s2, e2 = int(pos_matches[1].group("s")), int(pos_matches[1].group("e"))
            matches = int(h.group("matches"))
            total   = int(h.group("total"))
            score   = int(h.group("score"))
            gaps    = int(h.group("gaps"))
            seq1 = pos_matches[0].group("seq").replace(" ", "").replace("-", "")
            seq2 = pos_matches[1].group("seq").replace(" ", "").replace("-", "")
        except (ValueError, AttributeError):
            i = max(j, i + 1)
            continue
        left_start  = min(s1, e1);  left_end   = max(s1, e1)
        right_start = min(s2, e2);  right_end  = max(s2, e2)
        tir_length  = left_end - left_start + 1
        insert_size = max(0, right_start - left_end - 1)
        plv_span    = right_end - left_start + 1
        identity    = 100.0 * matches / total if total else 0.0
        pairs.append(TirPair(
            left_start=left_start, left_end=left_end,
            right_start=right_start, right_end=right_end,
            tir_length=tir_length, insert_size=insert_size, plv_span=plv_span,
            tir_identity=identity, score=score,
            matches=matches, total=total, gaps=gaps,
            left_seq=seq1, right_seq=seq2,
        ))
        i = j
    return pairs


def run_einverted_on_region(
    region_fa_path: Path, inv_out_path: Path, seq_out_path: Path, cfg: dict
) -> None:
    cmd = [
        "einverted",
        "-sequence",  str(region_fa_path),
        "-outfile",   str(inv_out_path),
        "-outseq",    str(seq_out_path),
        "-gap",       str(cfg["einverted_gap"]),
        "-threshold", str(cfg["einverted_threshold"]),
        "-match",     str(cfg["einverted_match"]),
        "-mismatch",  str(cfg["einverted_mismatch"]),
        "-maxrepeat", str(cfg["einverted_maxrepeat"]),
        "-auto",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(
            proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr
        )


def select_best_tir(
    pairs: List[TirPair],
    region_offset: int,
    cluster_start: int,
    cluster_end: int,
    mcp_intervals: List[Tuple[int, int]],
    cfg: dict,
) -> Tuple[Optional[TirPair], dict]:
    del cluster_start, cluster_end   # retained on signature for compatibility
    diag = dict(
        n_raw=len(pairs), n_pass_insert_size=0, n_pass_tir_length=0,
        n_pass_identity=0, n_pass_mcp_bracket=0, best_near_miss="",
    )
    require_mcp_bracket = bool(mcp_intervals)
    valid: List[TirPair] = []
    near_miss_score = -1.0
    for t in pairs:
        gleft_start  = t.left_start  + region_offset - 1
        gleft_end    = t.left_end    + region_offset - 1
        gright_start = t.right_start + region_offset - 1
        gright_end   = t.right_end   + region_offset - 1

        ok_insert   = cfg["tir_min_insert"] <= t.insert_size <= cfg["tir_max_insert"]
        ok_length   = cfg["tir_min_len"]    <= t.tir_length  <= cfg["tir_max_len"]
        ok_identity = t.tir_identity >= cfg["tir_min_id"]
        ok_mcp      = (
            any(gleft_start <= ms and me <= gright_end for ms, me in mcp_intervals)
            if require_mcp_bracket else True
        )

        if ok_insert:                                 diag["n_pass_insert_size"] += 1
        if ok_insert and ok_length:                   diag["n_pass_tir_length"] += 1
        if ok_insert and ok_length and ok_identity:   diag["n_pass_identity"] += 1
        if ok_insert and ok_length and ok_identity and ok_mcp:
            diag["n_pass_mcp_bracket"] += 1

        n_passed = sum([ok_insert, ok_length, ok_identity, ok_mcp])
        score = n_passed * 1e6 + t.tir_identity * t.tir_length
        if score > near_miss_score:
            near_miss_score = score
            diag["best_near_miss"] = (
                f"insert={t.insert_size:,}bp tir_len={t.tir_length}bp "
                f"id={t.tir_identity:.1f}% "
                f"left=[{gleft_start:,}-{gleft_end:,}] "
                f"right=[{gright_start:,}-{gright_end:,}] "
                f"passed=[size:{ok_insert},len:{ok_length},"
                f"id:{ok_identity},mcp:{ok_mcp}]"
            )
        if not (ok_insert and ok_length and ok_identity and ok_mcp):
            continue
        valid.append(TirPair(
            left_start=gleft_start, left_end=gleft_end,
            right_start=gright_start, right_end=gright_end,
            tir_length=t.tir_length, insert_size=t.insert_size, plv_span=t.plv_span,
            tir_identity=t.tir_identity, score=t.score,
            matches=t.matches, total=t.total, gaps=t.gaps,
            left_seq=t.left_seq, right_seq=t.right_seq,
        ))

    if not valid:
        return None, diag
    valid.sort(key=lambda x: (x.tir_identity, x.tir_length, x.insert_size), reverse=True)
    return valid[0], diag


# ── Stage 4: TSD detection ────────────────────────────────────────────────────
def find_tsd(left_flank: str, right_flank: str, k_min: int, k_max: int, max_slide: int) -> Optional[Tsd]:
    left, right = left_flank.upper(), right_flank.upper()
    best: Optional[Tsd] = None
    for k in range(k_max, k_min - 1, -1):
        if k > len(left) or k > len(right):
            continue
        max_mm = 0 if k <= 5 else (1 if k <= 8 else 2)
        for sl in range(0, max_slide + 1):
            if k + sl > len(left):
                break
            for sr in range(0, max_slide + 1):
                if k + sr > len(right):
                    break
                left_kmer  = left[len(left) - k - sl: len(left) - sl] if sl > 0 else left[-k:]
                right_kmer = right[sr: sr + k]
                if "N" in left_kmer or "N" in right_kmer:
                    continue
                mm = sum(1 for a, b in zip(left_kmer, right_kmer) if a != b)
                if mm <= max_mm:
                    identity = 100.0 * (k - mm) / k
                    cand = Tsd(
                        sequence_left=left_kmer, sequence_right=right_kmer,
                        length=k, mismatches=mm, identity=identity,
                        left_shift=sl, right_shift=sr,
                    )
                    if best is None or (cand.length, cand.identity) > (best.length, best.identity):
                        best = cand
        if best is not None and best.length == k:
            return best
    return best


# ── Stage 5: GC content ───────────────────────────────────────────────────────
def gc_in_windows(seq: str, window: int) -> np.ndarray:
    if not seq:
        return np.zeros(0, dtype=np.float32)
    arr = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    is_gc    = ((arr == ord("G")) | (arr == ord("C"))).astype(np.float32)
    is_valid = is_gc + ((arr == ord("A")) | (arr == ord("T"))).astype(np.float32)
    n = len(arr) // window
    if n == 0:
        return np.zeros(0, dtype=np.float32)
    gc_sum    = is_gc[:n*window].reshape(n, window).sum(axis=1)
    valid_sum = is_valid[:n*window].reshape(n, window).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(valid_sum > 0, 100.0 * gc_sum / valid_sum, np.nan)
    return pct.astype(np.float32)


def gc_mean_of_seq(seq: str) -> float:
    if not seq:
        return float("nan")
    arr  = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    is_gc = ((arr == ord("G")) | (arr == ord("C"))).sum()
    is_at = ((arr == ord("A")) | (arr == ord("T"))).sum()
    valid = is_gc + is_at
    return float(100.0 * is_gc / valid) if valid > 0 else float("nan")


# ── Stage 5b: TIR-less span control (anchor-driven + GC refinement) ──────────
def _orf_coding_bp(orfs: List[Orf], region_start: int, region_end: int) -> int:
    """Return the total ORF base pairs (annotated or unannotated) inside
    [region_start, region_end] (1-based inclusive). Overlapping ORFs are
    summed without merging, which is fine here as a relative density proxy.
    """
    if region_end < region_start:
        return 0
    total = 0
    for o in orfs:
        s = max(o.start, region_start)
        e = min(o.end,   region_end)
        if e >= s:
            total += e - s + 1
    return total


def trim_tirless_span(
    cluster_orfs: List[Orf],
    anchor_families: frozenset,
    edge_gap: int,
    max_length: int,
) -> Tuple[int, int, List[str]]:
    """Compute a robust span for a TIR-less PLV candidate.

    Three-phase, anchor-driven design.

    Phase 1 — Anchor bounding box.
        Span starts as [min anchor start, max anchor end] using ONLY ORFs
        whose family is in `anchor_families` (MCP, mCP, pPolB,
        A32-pATPase, RVE-INT). Supportive PLV families and unannotated
        ORFs cannot define the initial edges — this prevents host-prone
        domains (e.g. Y-rec, GIY-YIG) from over-extending the span.

    Phase 2 — Outward extension by adjacency.
        Walk outward from each anchor-box edge through the cluster ORFs
        sorted by position. An ORF (annotated OR unannotated) is pulled
        in if its gap to the current edge is <= `edge_gap`. The walk
        stops at the first ORF whose gap exceeds the threshold. This
        captures supportive PLV ORFs and viral-neighbourhood unannotated
        ORFs while excluding distant host genes.

    Phase 3 — Length cap that preserves all anchors.
        If the span exceeds `max_length`, anchors are NEVER evicted.
        Slack (max_length − anchor_span_length) is allocated between the
        two flanks proportional to local coding density (sum of ORF bp
        per flank length). The more PLV-like flank therefore keeps more
        of its content. If the anchor span itself already exceeds
        max_length, the span is forced to the anchor box and a warning
        is emitted.

    Parameters
    ----------
    cluster_orfs : list of Orf
        ALL ORFs (annotated + unannotated) within the original cluster
        bounds. Phase 2 considers any of these for outward extension.
    anchor_families : frozenset
        The set of family labels treated as anchors (typically
        ANCHOR_FAMILIES_BASE plus mCP if its HMM was supplied).
    edge_gap : int
        Maximum allowed gap (bp) for Phase 2 outward extension.
    max_length : int
        Hard cap on the final span length (bp). May be exceeded only when
        the anchor span itself already does (a warning is emitted).

    Returns
    -------
    (span_start, span_end, notes) : (int, int, list of str)
        Genome-relative 1-based inclusive coordinates of the trimmed
        span, plus human-readable notes for the run log.
    """
    notes: List[str] = []

    if not cluster_orfs:
        return 0, 0, ["no ORFs in cluster; span undefined"]

    anchors = [o for o in cluster_orfs if o.family in anchor_families]
    if not anchors:
        # Should not normally occur (seeding requires MCP), but stay safe.
        annotated = [o for o in cluster_orfs if o.family is not None]
        if not annotated:
            return 0, 0, ["no annotated/anchor ORFs; span undefined"]
        ann_sorted = sorted(annotated, key=lambda o: o.start)
        return (
            ann_sorted[0].start,
            ann_sorted[-1].end,
            ["no anchor ORFs found; falling back to annotated bounding box"],
        )

    anchors_sorted = sorted(anchors, key=lambda o: o.start)
    anchor_min = anchors_sorted[0].start
    anchor_max = anchors_sorted[-1].end
    span_start, span_end = anchor_min, anchor_max
    anchor_fams_str = ",".join(sorted({o.family for o in anchors_sorted}))
    notes.append(
        f"Phase 1: anchor bounding box {span_start:,}-{span_end:,} bp "
        f"({len(anchors_sorted)} anchor ORF(s); families: {anchor_fams_str})"
    )

    # ── Phase 2: outward extension by adjacency ──────────────────────────────
    cluster_sorted = sorted(cluster_orfs, key=lambda o: o.start)

    # Leftward (iterate ORFs strictly left of current edge, in reverse).
    left_candidates = [o for o in cluster_sorted if o.end < span_start]
    n_left = 0
    for o in reversed(left_candidates):
        gap = span_start - o.end - 1
        if gap > edge_gap:
            break
        tag = o.family if o.family else "unannotated"
        notes.append(
            f"Phase 2: extended left to include {o.orf_id} ({tag}; "
            f"gap={gap:,} bp <= {edge_gap:,} bp)"
        )
        span_start = o.start
        n_left += 1

    # Rightward (iterate ORFs strictly right of current edge).
    right_candidates = [o for o in cluster_sorted if o.start > span_end]
    n_right = 0
    for o in right_candidates:
        gap = o.start - span_end - 1
        if gap > edge_gap:
            break
        tag = o.family if o.family else "unannotated"
        notes.append(
            f"Phase 2: extended right to include {o.orf_id} ({tag}; "
            f"gap={gap:,} bp <= {edge_gap:,} bp)"
        )
        span_end = o.end
        n_right += 1

    if n_left == 0 and n_right == 0:
        notes.append("Phase 2: no qualifying flanking ORFs; span unchanged")

    # ── Phase 3: length cap, preserving every anchor ─────────────────────────
    span_len = span_end - span_start + 1
    if span_len <= max_length:
        return span_start, span_end, notes

    anchor_span_len = anchor_max - anchor_min + 1
    if anchor_span_len > max_length:
        notes.append(
            f"Phase 3: WARNING — anchor span ({anchor_span_len:,} bp) exceeds "
            f"max_length ({max_length:,} bp); preserving all anchors regardless"
        )
        return anchor_min, anchor_max, notes

    slack          = max_length - anchor_span_len
    left_flank_len  = anchor_min - span_start
    right_flank_len = span_end   - anchor_max

    left_density_bp  = _orf_coding_bp(cluster_orfs, span_start,    anchor_min - 1)
    right_density_bp = _orf_coding_bp(cluster_orfs, anchor_max + 1, span_end)
    left_density  = left_density_bp  / left_flank_len  if left_flank_len  > 0 else 0.0
    right_density = right_density_bp / right_flank_len if right_flank_len > 0 else 0.0

    total_density = left_density + right_density
    if total_density > 0:
        left_alloc  = int(slack * left_density / total_density)
        right_alloc = slack - left_alloc
    else:
        left_alloc  = slack // 2
        right_alloc = slack - left_alloc

    # Cap allocations by actual flank size, then redistribute leftover.
    left_alloc  = min(left_alloc,  left_flank_len)
    right_alloc = min(right_alloc, right_flank_len)
    leftover = slack - left_alloc - right_alloc
    if leftover > 0 and left_alloc < left_flank_len:
        extra = min(leftover, left_flank_len - left_alloc)
        left_alloc += extra
        leftover -= extra
    if leftover > 0 and right_alloc < right_flank_len:
        extra = min(leftover, right_flank_len - right_alloc)
        right_alloc += extra

    new_start = anchor_min - left_alloc
    new_end   = anchor_max + right_alloc
    notes.append(
        f"Phase 3: capped {span_len:,} bp -> {new_end - new_start + 1:,} bp "
        f"[{new_start:,}-{new_end:,}]; flank slack allocated by coding density "
        f"(left={left_alloc:,} bp at {left_density:.3f} bp/bp, "
        f"right={right_alloc:,} bp at {right_density:.3f} bp/bp); "
        f"all {len(anchors_sorted)} anchor(s) preserved"
    )
    return new_start, new_end, notes


def refine_span_by_gc(
    span_start: int,
    span_end: int,
    anchor_min: int,
    anchor_max: int,
    seq: str,
    seq_offset: int,
    contig_length: int,
    window: int,
    flank_size: int,
    min_delta_pct: float,
) -> Tuple[int, int, List[str]]:
    """Snap span boundaries inward at GC% discontinuities.

    Compares mean GC% in the anchor span (likely PLV-composition) with
    mean GC% in the immediate flank just outside the current span (likely
    host-composition). When the absolute difference exceeds
    `min_delta_pct`, walks inward window-by-window from each edge and
    trims contiguous host-like windows (i.e. windows whose GC% is closer
    to the flank GC than to the anchor GC). Never trims past an anchor.

    Conservative by design: only shrinks the span, never extends it. If
    no significant compositional discontinuity is detected, returns the
    input boundaries unchanged.

    Parameters
    ----------
    span_start, span_end : int
        Current span (1-based inclusive, genome coordinates).
    anchor_min, anchor_max : int
        Anchor bounding box; refinement never crosses these.
    seq : str
        Substring of the contig sequence covering at least
        [span_start − flank_size, span_end + flank_size] when possible.
    seq_offset : int
        Genome-relative 1-based start position of `seq`.
    contig_length : int
        Length of the parent contig (used only for clamping; flanks at
        contig ends just yield NaN GC and are skipped gracefully).
    window : int
        Sliding-window size (bp) for GC% computation.
    flank_size : int
        Bases of immediate flank used to estimate host GC.
    min_delta_pct : float
        Minimum |anchor GC − flank GC| (percentage points) required to
        trigger refinement on that side.
    """
    notes: List[str] = []

    def _slice(start: int, end: int) -> str:
        # Genome 1-based inclusive -> local seq indices.
        local_s = max(0, start - seq_offset)
        local_e = min(len(seq), end - seq_offset + 1)
        if local_e <= local_s:
            return ""
        return seq[local_s:local_e]

    def _gc_pct(s: str) -> float:
        if not s:
            return float("nan")
        u = s.upper()
        gc = u.count("G") + u.count("C")
        at = u.count("A") + u.count("T")
        n  = gc + at
        return 100.0 * gc / n if n > 0 else float("nan")

    gc_anchor = _gc_pct(_slice(anchor_min, anchor_max))
    if np.isnan(gc_anchor):
        return span_start, span_end, ["GC refinement skipped: anchor GC undefined"]

    gc_left_flank  = _gc_pct(_slice(span_start - flank_size, span_start - 1))
    gc_right_flank = _gc_pct(_slice(span_end + 1, span_end + flank_size))

    new_start, new_end = span_start, span_end

    # ── Left edge ────────────────────────────────────────────────────────────
    if not np.isnan(gc_left_flank):
        delta_left = abs(gc_anchor - gc_left_flank)
        if delta_left >= min_delta_pct:
            pos = span_start
            while pos + window - 1 < anchor_min:
                gc_w = _gc_pct(_slice(pos, pos + window - 1))
                if np.isnan(gc_w):
                    break
                d_anchor = abs(gc_w - gc_anchor)
                d_flank  = abs(gc_w - gc_left_flank)
                if d_flank < d_anchor:
                    pos += window
                else:
                    break
            if pos > span_start:
                notes.append(
                    f"GC refinement: trimmed left edge {span_start:,} -> {pos:,} "
                    f"(anchor GC={gc_anchor:.1f}%, "
                    f"left-flank GC={gc_left_flank:.1f}%, "
                    f"delta={delta_left:.1f} pp)"
                )
                new_start = pos

    # ── Right edge ───────────────────────────────────────────────────────────
    if not np.isnan(gc_right_flank):
        delta_right = abs(gc_anchor - gc_right_flank)
        if delta_right >= min_delta_pct:
            pos = span_end
            while pos - window + 1 > anchor_max:
                gc_w = _gc_pct(_slice(pos - window + 1, pos))
                if np.isnan(gc_w):
                    break
                d_anchor = abs(gc_w - gc_anchor)
                d_flank  = abs(gc_w - gc_right_flank)
                if d_flank < d_anchor:
                    pos -= window
                else:
                    break
            if pos < span_end:
                notes.append(
                    f"GC refinement: trimmed right edge {span_end:,} -> {pos:,} "
                    f"(anchor GC={gc_anchor:.1f}%, "
                    f"right-flank GC={gc_right_flank:.1f}%, "
                    f"delta={delta_right:.1f} pp)"
                )
                new_end = pos

    if not notes:
        notes.append("GC refinement: no significant discontinuity at edges")
    return new_start, new_end, notes


# ── Stage 6 helper: output writers ───────────────────────────────────────────
def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def stop_codon_of_orf(orf: Orf, fa: "pyfastx.Fasta", contig_length: int) -> Optional[str]:
    if orf.partial:
        return None
    if orf.strand >= 0:
        if orf.end < 3 or orf.end > contig_length:
            return None
        seq = fa.fetch(orf.contig, (orf.end - 2, orf.end))
    else:
        if orf.start < 1 or orf.start + 2 > contig_length:
            return None
        seq = reverse_complement(fa.fetch(orf.contig, (orf.start, orf.start + 2)))
    seq = seq.upper()
    return seq if len(seq) == 3 else None


def compute_plv_stats(plv: Plv, fa: "pyfastx.Fasta") -> dict:
    n_orfs = len(plv.orfs)
    coding_bases = 0
    if plv.length > 0 and plv.orfs:
        intervals = [(max(o.start, plv.start), min(o.end, plv.end)) for o in plv.orfs]
        intervals = [(s, e) for s, e in intervals if e >= s]
        intervals.sort()
        merged_ivs: List[Tuple[int, int]] = []
        for s, e in intervals:
            if merged_ivs and s <= merged_ivs[-1][1] + 1:
                merged_ivs[-1] = (merged_ivs[-1][0], max(merged_ivs[-1][1], e))
            else:
                merged_ivs.append((s, e))
        coding_bases = sum(e - s + 1 for s, e in merged_ivs)
    coding_pct    = min(100.0 * coding_bases / plv.length, 100.0) if plv.length > 0 else 0.0
    noncoding_pct = 100.0 - coding_pct

    stop_counts = {"TAG": 0, "TAA": 0, "TGA": 0, "OTHER": 0}
    n_complete = 0
    for o in plv.orfs:
        sc = stop_codon_of_orf(o, fa, plv.contig_length)
        if sc is None:
            continue
        n_complete += 1
        stop_counts[sc if sc in stop_counts else "OTHER"] += 1
    stop_pcts = (
        {k: round(100.0 * v / n_complete, 2) for k, v in stop_counts.items()}
        if n_complete > 0 else {k: 0.0 for k in stop_counts}
    )
    return dict(
        n_orfs=n_orfs, n_orfs_complete=n_complete,
        coding_density_pct=round(coding_pct, 3),
        noncoding_density_pct=round(noncoding_pct, 3),
        stop_tag_pct=stop_pcts["TAG"], stop_taa_pct=stop_pcts["TAA"],
        stop_tga_pct=stop_pcts["TGA"], stop_other_pct=stop_pcts["OTHER"],
    )


def format_tsd(tsd: Optional[Tsd]) -> Tuple[str, str]:
    if tsd is None:
        return "NA", "NODETECT"
    details = (
        f"[LEN={tsd.length},D1={tsd.left_shift},D2={tsd.right_shift},"
        f"SEQ={tsd.sequence_left}]"
    )
    return details, ("PERFECT" if tsd.mismatches == 0 else "IMPERFECT")


def write_run_summary_txt(
    plvs: List[Plv], path: Path, prefix: str, genome_path: Path, elapsed_seconds: float
) -> None:
    def fmt_pct(n, d):
        return f"{100.0 * n / d:.1f}%" if d > 0 else "n/a"

    def fmt_elapsed(s):
        if s >= 3600: return f"{s / 3600:.2f} h"
        if s >= 60:   return f"{s / 60:.2f} min"
        return f"{s:.1f} s"

    n_total    = len(plvs)
    n_with_tir = sum(1 for p in plvs if p.has_tir)
    lines: List[str] = ["Summary Statistics", "-" * 60,
        f"Input genome:      {genome_path}",
        f"Elapsed:           {fmt_elapsed(elapsed_seconds)}",
        f"Candidate PLVs:    {n_total}",
    ]
    if n_total == 0:
        lines.append("No PLV candidates were reported. See run.log for details.")
        path.write_text("\n".join(lines) + "\n")
        return

    n_no_tir = n_total - n_with_tir
    lines += [
        f"  With TIR:        {n_with_tir} ({fmt_pct(n_with_tir, n_total)})",
        f"  Without TIR:     {n_no_tir} ({fmt_pct(n_no_tir, n_total)})",
    ]

    tir_plvs = [p for p in plvs if p.has_tir]
    if tir_plvs:
        n_perfect   = sum(1 for p in tir_plvs if p.tsd is not None and p.tsd.mismatches == 0)
        n_imperfect = sum(1 for p in tir_plvs if p.tsd is not None and p.tsd.mismatches > 0)
        n_nodetect  = sum(1 for p in tir_plvs if p.tsd is None)
        lines += [
            "TSD detection (TIR-bearing PLVs):",
            f"  PERFECT:         {n_perfect} ({fmt_pct(n_perfect, len(tir_plvs))})",
            f"  IMPERFECT:       {n_imperfect} ({fmt_pct(n_imperfect, len(tir_plvs))})",
            f"  NODETECT:        {n_nodetect} ({fmt_pct(n_nodetect, len(tir_plvs))})",
        ]

    lengths = np.array([p.length for p in plvs])
    lines += [
        "PLV length (bp):",
        f"  Min:             {int(lengths.min()):,}",
        f"  Max:             {int(lengths.max()):,}",
        f"  Mean:            {int(lengths.mean()):,}",
        f"  Median:          {int(np.median(lengths)):,}",
    ]

    fam_counter: Dict[str, int] = {}
    for p in plvs:
        for f in p.families_present:
            fam_counter[f] = fam_counter.get(f, 0) + 1
    top_fams = sorted(fam_counter.items(), key=lambda kv: kv[1], reverse=True)[:6]
    top_str  = ", ".join(f"{n}({c})" for n, c in top_fams)
    n_fams_per_plv = [p.n_families for p in plvs]
    lines += [
        "Marker families:",
        f"  Mean per PLV:    {np.mean(n_fams_per_plv):.1f}",
        f"  Most frequent:   {top_str}",
    ]

    if tir_plvs:
        tl = np.array([p.tir.tir_length for p in tir_plvs])
        lines += [
            "TIR length (bp):",
            f"  Min:             {int(tl.min()):,}",
            f"  Max:             {int(tl.max()):,}",
            f"  Mean:            {int(tl.mean()):,}",
            f"  Median:          {int(np.median(tl)):,}",
        ]

    gc_vals = np.array([p.gc_plv for p in plvs if not np.isnan(p.gc_plv)])
    if len(gc_vals) > 0:
        lines += ["GC content (%):", f"  Mean PLV:        {gc_vals.mean():.1f}"]

    lines.append(f"Contigs carrying PLVs: {len({p.contig for p in plvs})}")
    path.write_text("\n".join(lines) + "\n")


def write_plv_tsv(
    plvs: List[Plv], fa: Optional["pyfastx.Fasta"], path: Path, prefix: str
) -> None:
    rows = []
    for p in plvs:
        orf_annots = []
        for o in p.orfs:
            fam = o.family if o.family else "Unannotated"
            low_conf = ""
            if o.family is None and o.best_pfam_label:
                low_conf = (
                    f"[best_hit:{o.best_pfam_label}|{o.best_pfam_acc}"
                    f"|bs={o.best_pfam_bitscore:.1f}|e={o.best_pfam_evalue:.2e}]"
                )
            orf_annots.append(f"{o.orf_id}={fam}{low_conf}")

        stats = compute_plv_stats(p, fa) if fa is not None else dict(
            n_orfs=len(p.orfs), n_orfs_complete=0,
            coding_density_pct=0.0, noncoding_density_pct=0.0,
            stop_tag_pct=0.0, stop_taa_pct=0.0,
            stop_tga_pct=0.0, stop_other_pct=0.0,
        )

        if p.tir is not None:
            tla, tle = p.tir.left_start,   p.tir.left_end
            tra, tre = p.tir.right_start,  p.tir.right_end
            tlr  = p.tir.left_start   - p.start + 1
            tler = p.tir.left_end     - p.start + 1
            trr  = p.tir.right_start  - p.start + 1
            trer = p.tir.right_end    - p.start + 1
        else:
            tla = tle = tra = tre = "NA"
            tlr = tler = trr = trer = "NA"

        tsd_details, tsd_conservation = format_tsd(p.tsd)

        rows.append(dict(
            plv_id=p.plv_id, prefix=prefix,
            contig=p.contig, contig_length=p.contig_length,
            start=p.start, end=p.end, length=p.length,
            has_tir=("yes" if p.has_tir else "no"),
            n_families=p.n_families, families=";".join(p.families_present),
            mcp_best_bitscore=round(p.mcp_best_bitscore, 2),
            n_orfs=stats["n_orfs"], n_orfs_complete=stats["n_orfs_complete"],
            coding_density_pct=stats["coding_density_pct"],
            noncoding_density_pct=stats["noncoding_density_pct"],
            stop_tag_pct=stats["stop_tag_pct"], stop_taa_pct=stats["stop_taa_pct"],
            stop_tga_pct=stats["stop_tga_pct"], stop_other_pct=stats["stop_other_pct"],
            gc_plv_pct=(round(p.gc_plv, 2) if not np.isnan(p.gc_plv) else "NA"),
            tir_left_start_abs=tla,   tir_left_end_abs=tle,
            tir_right_start_abs=tra,  tir_right_end_abs=tre,
            tir_left_start_rel=tlr,   tir_left_end_rel=tler,
            tir_right_start_rel=trr,  tir_right_end_rel=trer,
            tir_left_length=(tler - tlr + 1 if p.tir else "NA"),
            tir_right_length=(trer - trr + 1 if p.tir else "NA"),
            tir_score=(p.tir.score if p.tir else "NA"),
            tir_matches=(p.tir.matches if p.tir else "NA"),
            tir_total=(p.tir.total if p.tir else "NA"),
            tir_identity_pct=(round(p.tir.tir_identity, 2) if p.tir else "NA"),
            tir_gaps=(p.tir.gaps if p.tir else "NA"),
            tir_orient_left=("+" if p.tir else "NA"),
            tir_orient_right=("-" if p.tir else "NA"),
            insert_length=(p.tir.insert_size if p.tir else "NA"),
            tsd_length=(p.tsd.length if p.tsd else 0),
            tsd_identity_pct=(round(p.tsd.identity, 2) if p.tsd else 0.0),
            tsd_mismatches=(p.tsd.mismatches if p.tsd else 0),
            tsd_left_seq=(p.tsd.sequence_left if p.tsd else ""),
            tsd_right_seq=(p.tsd.sequence_right if p.tsd else ""),
            tsd_left_dist=(p.tsd.left_shift if p.tsd else 0),
            tsd_right_dist=(p.tsd.right_shift if p.tsd else 0),
            tsd_details=tsd_details, tsd_conservation=tsd_conservation,
            orf_annotations="|".join(orf_annots),
        ))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def write_multifasta(plvs: List[Plv], fa_index: pyfastx.Fasta, path: Path, prefix: str) -> None:
    with open(path, "w") as fh:
        for p in plvs:
            seq     = fa_index.fetch(p.contig, (p.start, p.end))
            tsd_str = (p.tsd.sequence_left if p.tsd else "")
            tir_str = "yes" if p.has_tir else "no"
            header  = (
                f">{p.plv_id} prefix={prefix} contig={p.contig} "
                f"start={p.start} end={p.end} length={p.length} has_tir={tir_str} "
            )
            if p.tir:
                header += (
                    f"tirL={p.tir.left_start}..{p.tir.left_end} "
                    f"tirR={p.tir.right_start}..{p.tir.right_end} "
                )
            header += f"tsd={tsd_str or 'NA'} gc_plv_pct={p.gc_plv:.2f}%"
            fh.write(header + "\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def write_gff3(plvs: List[Plv], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for p in plvs:
            tir_note = "has_tir=yes" if p.has_tir else "has_tir=no"
            attrs = (
                f"ID={p.plv_id};Name={p.plv_id};"
                f"length={p.length};families={','.join(p.families_present)};"
                f"{tir_note}"
            )
            fh.write(
                f"{p.contig}\tfindPLV\tmobile_genetic_element\t"
                f"{p.start}\t{p.end}\t.\t+\t.\t{attrs}\n"
            )
            if p.tir is not None:
                fh.write(
                    f"{p.contig}\tfindPLV\tterminal_inverted_repeat\t"
                    f"{p.tir.left_start}\t{p.tir.left_end}\t{p.tir.score}\t+\t.\t"
                    f"ID={p.plv_id}_tirL;Parent={p.plv_id};"
                    f"identity={p.tir.tir_identity:.1f}\n"
                )
                fh.write(
                    f"{p.contig}\tfindPLV\tterminal_inverted_repeat\t"
                    f"{p.tir.right_start}\t{p.tir.right_end}\t{p.tir.score}\t-\t.\t"
                    f"ID={p.plv_id}_tirR;Parent={p.plv_id};"
                    f"identity={p.tir.tir_identity:.1f}\n"
                )
            if p.tsd is not None and p.tir is not None:
                ltsd_end   = p.tir.left_start  - 1 - p.tsd.left_shift
                ltsd_start = ltsd_end - p.tsd.length + 1
                rtsd_start = p.tir.right_end   + 1 + p.tsd.right_shift
                rtsd_end   = rtsd_start + p.tsd.length - 1
                if ltsd_start >= 1:
                    fh.write(
                        f"{p.contig}\tfindPLV\ttarget_site_duplication\t"
                        f"{ltsd_start}\t{ltsd_end}\t.\t+\t.\t"
                        f"ID={p.plv_id}_tsdL;Parent={p.plv_id};"
                        f"sequence={p.tsd.sequence_left};length={p.tsd.length};"
                        f"mismatches={p.tsd.mismatches}\n"
                    )
                if rtsd_end <= p.contig_length:
                    fh.write(
                        f"{p.contig}\tfindPLV\ttarget_site_duplication\t"
                        f"{rtsd_start}\t{rtsd_end}\t.\t+\t.\t"
                        f"ID={p.plv_id}_tsdR;Parent={p.plv_id};"
                        f"sequence={p.tsd.sequence_right};length={p.tsd.length};"
                        f"mismatches={p.tsd.mismatches}\n"
                    )
            for o in p.orfs:
                family    = o.family if o.family else "Unannotated"
                src       = o.family_source if o.family_source else "."
                strand_ch = "+" if o.strand >= 0 else "-"
                extra = ""
                if o.family is None and o.best_pfam_label:
                    extra = (
                        f";best_pfam_hit={o.best_pfam_label}"
                        f";best_pfam_acc={o.best_pfam_acc}"
                        f";best_pfam_bitscore={o.best_pfam_bitscore:.1f}"
                        f";best_pfam_evalue={o.best_pfam_evalue:.2e}"
                    )
                fh.write(
                    f"{p.contig}\tfindPLV\tCDS\t{o.start}\t{o.end}\t"
                    f"{o.family_bitscore:.1f}\t{strand_ch}\t0\t"
                    f"ID={o.orf_id};Parent={p.plv_id};family={family};"
                    f"source={src};evalue={o.family_evalue:.2e}{extra}\n"
                )


def write_gggenomes(plvs: List[Plv], seqs_path: Path, genes_path: Path) -> None:
    seqs_rows = [dict(
        seq_id=p.plv_id, length=p.length, contig=p.contig,
        contig_start=p.start, contig_end=p.end,
        has_tir=("yes" if p.has_tir else "no"),
        n_families=p.n_families,
    ) for p in plvs]
    pd.DataFrame(seqs_rows).to_csv(seqs_path, sep="\t", index=False)

    gene_rows = []
    for p in plvs:
        if p.tir is not None:
            gene_rows.append(dict(
                seq_id=p.plv_id, start=1, end=p.tir.tir_length,
                strand="+", feat_id=f"{p.plv_id}_tirL", type="TIR", name="TIR_L",
            ))
            r_local_s = p.tir.right_start - p.start + 1
            r_local_e = p.tir.right_end   - p.start + 1
            gene_rows.append(dict(
                seq_id=p.plv_id, start=r_local_s, end=r_local_e,
                strand="-", feat_id=f"{p.plv_id}_tirR", type="TIR", name="TIR_R",
            ))
        for o in p.orfs:
            local_start = max(1,         o.start - p.start + 1)
            local_end   = min(p.length,  o.end   - p.start + 1)
            if o.family:
                family = o.family
            elif o.best_pfam_label:
                family = f"{o.best_pfam_label}(low_conf)"
            else:
                family = "Unannotated"
            gene_rows.append(dict(
                seq_id=p.plv_id, start=local_start, end=local_end,
                strand=("+" if o.strand >= 0 else "-"),
                feat_id=o.orf_id, type="CDS", name=family,
            ))
    pd.DataFrame(gene_rows).to_csv(genes_path, sep="\t", index=False)


def write_bedgraph(
    plvs: List[Plv], fa_index: pyfastx.Fasta, contig_lengths: Dict[str, int],
    path: Path, gc_window: int, flank: int,
) -> None:
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for p in plvs:
        s = max(1, p.start - flank)
        e = min(contig_lengths[p.contig], p.end + flank)
        by_contig.setdefault(p.contig, []).append((s, e))
    for k in by_contig:
        ivs = sorted(by_contig[k])
        merged = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        by_contig[k] = merged

    with open(path, "w") as fh:
        fh.write(
            f"track type=bedGraph name=GC_pct "
            f"description=\"GC% in {gc_window} bp windows\"\n"
        )
        for contig, regions in by_contig.items():
            for rstart, rend in regions:
                seq = fa_index.fetch(contig, (rstart, rend))
                pct = gc_in_windows(seq, gc_window)
                for i, val in enumerate(pct):
                    if np.isnan(val):
                        continue
                    wstart = rstart + i * gc_window - 1
                    wend   = wstart + gc_window
                    fh.write(f"{contig}\t{wstart}\t{wend}\t{val:.2f}\n")


def build_family_table() -> Dict[str, Tuple[str, bool]]:
    return {acc: (label, seed) for acc, label, seed in DEFAULT_PLV_FAMILIES}


# ── PLV deduplication ─────────────────────────────────────────────────────────
def deduplicate_plvs(
    plvs: List[Plv],
    min_reciprocal_overlap: float,
) -> List[Plv]:
    """Remove redundant overlapping PLV calls on the same contig.

    When two PLVs share > `min_reciprocal_overlap` of the shorter element's
    length, the weaker call (lower MCP bitscore; families as tiebreak) is
    discarded.  The ranking ensures we always keep the best-supported call.
    """
    if len(plvs) <= 1:
        return plvs

    # Sort: best-supported first (highest MCP bitscore, then most families).
    ranked = sorted(
        plvs,
        key=lambda p: (p.mcp_best_bitscore, p.n_families),
        reverse=True,
    )
    keep = [True] * len(ranked)
    for i in range(len(ranked)):
        if not keep[i]:
            continue
        for j in range(i + 1, len(ranked)):
            if not keep[j]:
                continue
            a, b = ranked[i], ranked[j]
            if a.contig != b.contig:
                continue
            overlap = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
            if overlap <= 0:
                continue
            shorter = min(a.length, b.length)
            if shorter > 0 and overlap / shorter >= min_reciprocal_overlap:
                keep[j] = False
                log_warning(
                    f"Deduplication: dropping {b.contig}:{b.start:,}-{b.end:,} "
                    f"({b.length:,} bp, bitscore={b.mcp_best_bitscore:.1f}) — "
                    f"{overlap / shorter * 100:.0f}% reciprocal overlap with "
                    f"{a.contig}:{a.start:,}-{a.end:,} (kept, "
                    f"bitscore={a.mcp_best_bitscore:.1f})"
                )
    kept_plvs = [p for p, k in zip(ranked, keep) if k]
    n_dropped = len(plvs) - len(kept_plvs)
    if n_dropped:
        log_info(f"Deduplication: removed {n_dropped} redundant PLV call(s)")
    return kept_plvs


# ── Stage 3-5 parallel worker ─────────────────────────────────────────────────
def _process_cluster(task: dict) -> dict:
    """Process one marker cluster: einverted -> TIR -> TSD -> GC.

    TIR-bearing PLVs
    ----------------
    Span = TIR left_start .. TIR right_end (enforced by select_best_tir's
    MCP-bracketing filter). No span refinement is applied.

    TIR-less PLVs — anchor-driven span with GC refinement
    -----------------------------------------------------
    1. trim_tirless_span() runs the three-phase logic:
       Phase 1 = anchor bounding box;
       Phase 2 = outward extension while gap <= edge_gap_trim;
       Phase 3 = max-length cap that preserves all anchors and
       distributes flank slack by local coding density.
    2. refine_span_by_gc() then conservatively snaps the boundaries
       inward at GC% discontinuities (only when |anchor GC − flank GC|
       >= gc_min_delta_pct). Never trims past an anchor.

    All candidates (TIR or not) that survive the length and MCP-in-span
    filters are retained; the has_tir flag records which case applies.
    """
    ci     = task["cluster_index"]
    contig = task["contig"]
    cstart = task["cluster_start"]
    cend   = task["cluster_end"]
    clen   = task["contig_length"]
    cfg    = task["cfg"]
    genome_path  = task["genome_path"]
    contig_orfs  = task["contig_orfs"]

    fa = pyfastx.Fasta(str(genome_path), build_index=True, uppercase=True)

    rstart     = max(1, cstart - cfg["flank_for_tir"])
    rend       = min(clen, cend + cfg["flank_for_tir"])
    region_seq = fa.fetch(contig, (rstart, rend))

    mcp_intervals: List[Tuple[int, int]] = [
        (o.start, o.end)
        for o in contig_orfs
        if o.family == MCP_LABEL and o.start <= cend and o.end >= cstart
    ]

    log_msgs: List[Tuple[str, str]] = []

    # ── Stage 3: TIR detection ────────────────────────────────────────────────
    best_tir: Optional[TirPair] = None
    tir_note: str = ""

    with tempfile.TemporaryDirectory(prefix="findPLV_") as tmp:
        tmp = Path(tmp)
        region_fa = tmp / f"region_{ci:04d}.fa"
        inv_out   = tmp / f"region_{ci:04d}.inv"
        seq_out   = tmp / f"region_{ci:04d}.fas"
        with open(region_fa, "w") as fh:
            fh.write(f">region_{ci:04d}\n")
            for i in range(0, len(region_seq), 80):
                fh.write(region_seq[i:i + 80] + "\n")
        try:
            run_einverted_on_region(region_fa, inv_out, seq_out, cfg)
        except FileNotFoundError:
            return dict(
                status="fatal", cluster_index=ci, contig=contig,
                message="einverted not found in PATH inside worker process",
            )
        except subprocess.CalledProcessError as e:
            tir_note = f"einverted failed: {e.stderr.strip()[:120]}"
            log_msgs.append(("warning",
                f"Cluster {ci} ({contig}:{cstart:,}-{cend:,}): "
                f"{tir_note}; retaining as TIR-less candidate"))
        else:
            pairs = parse_einverted(inv_out)
            best_tir, tir_diag = select_best_tir(
                pairs, region_offset=rstart,
                cluster_start=cstart, cluster_end=cend,
                mcp_intervals=mcp_intervals, cfg=cfg,
            )
            if best_tir is None:
                if tir_diag["n_raw"] == 0:
                    tir_note = "einverted reported no inverted-repeat pairs in this region"
                else:
                    funnel = (
                        f"{tir_diag['n_raw']} raw -> "
                        f"{tir_diag['n_pass_insert_size']} pass insert-size -> "
                        f"{tir_diag['n_pass_tir_length']} pass TIR-length -> "
                        f"{tir_diag['n_pass_identity']} pass identity -> "
                        f"{tir_diag['n_pass_mcp_bracket']} bracket MCP"
                    )
                    tir_note = (
                        f"all TIR candidates filtered [{funnel}]"
                        + (f"; best near-miss: {tir_diag['best_near_miss']}"
                           if tir_diag["best_near_miss"] else "")
                    )
                log_msgs.append(("info",
                    f"Cluster {ci} ({contig}:{cstart:,}-{cend:,}): "
                    f"{tir_note}; retaining as TIR-less candidate"))

    # ── PLV span determination ────────────────────────────────────────────────
    if best_tir is not None:
        # TIR case: span is exactly the TIR-bracketed region.
        plv_start = best_tir.left_start
        plv_end   = best_tir.right_end
    else:
        # TIR-less: anchor-driven three-phase span + GC refinement.
        anchor_families = cfg["anchor_families"]
        cluster_orfs_in = [
            o for o in contig_orfs if o.start >= cstart and o.end <= cend
        ]

        if not cluster_orfs_in:
            # Defensive fallback (should not occur after seeding).
            plv_start, plv_end = cstart, cend
        else:
            plv_start, plv_end, trim_notes = trim_tirless_span(
                cluster_orfs_in,
                anchor_families=anchor_families,
                edge_gap=cfg["edge_gap_trim"],
                max_length=cfg["max_plv_length"],
            )
            for note in trim_notes:
                log_msgs.append(("info",
                    f"Cluster {ci} ({contig}:{cstart:,}-{cend:,}): {note}"))

            # GC-based refinement only when we actually have anchors.
            anchors_in = [
                o for o in cluster_orfs_in if o.family in anchor_families
            ]
            if anchors_in and plv_end > plv_start:
                anchor_min_pos = min(a.start for a in anchors_in)
                anchor_max_pos = max(a.end   for a in anchors_in)
                gc_seq_start = max(1, plv_start - cfg["gc_flank_search"] - 10)
                gc_seq_end   = min(clen, plv_end + cfg["gc_flank_search"] + 10)
                gc_seq = fa.fetch(contig, (gc_seq_start, gc_seq_end))
                plv_start, plv_end, gc_notes = refine_span_by_gc(
                    span_start=plv_start, span_end=plv_end,
                    anchor_min=anchor_min_pos, anchor_max=anchor_max_pos,
                    seq=gc_seq, seq_offset=gc_seq_start,
                    contig_length=clen,
                    window=cfg["gc_window"],
                    flank_size=cfg["gc_flank_search"],
                    min_delta_pct=cfg["gc_min_delta_pct"],
                )
                for note in gc_notes:
                    log_msgs.append(("info",
                        f"Cluster {ci} ({contig}:{cstart:,}-{cend:,}): {note}"))

    # ── Length filter ─────────────────────────────────────────────────────────
    plv_length = plv_end - plv_start + 1
    if plv_length < cfg["min_plv_length"]:
        return dict(
            status="skip", cluster_index=ci, contig=contig, log_msgs=log_msgs,
            message=(
                f"Cluster {ci} ({contig}:{cstart:,}-{cend:,}): PLV span "
                f"{plv_length:,} bp < {cfg['min_plv_length']:,} bp threshold; discarding"
            ),
        )

    # ── ORFs within the final PLV span ────────────────────────────────────────
    plv_orfs = [o for o in contig_orfs if o.start >= plv_start and o.end <= plv_end]
    plv_orfs.sort(key=lambda x: x.start)
    fams_in_plv = {o.family for o in plv_orfs if o.family is not None}

    # Safety net: MCP must be inside the reported span.
    if cfg["require_mcp"] and MCP_LABEL not in fams_in_plv:
        return dict(
            status="skip", cluster_index=ci, contig=contig, log_msgs=log_msgs,
            message=f"Cluster {ci}: MCP ORF outside final PLV boundaries; discarding",
        )

    mcp_best = max(
        (o.family_bitscore for o in plv_orfs if o.family == MCP_LABEL),
        default=0.0,
    )

    # ── Stage 4: TSD (only with known TIR boundaries) ────────────────────────
    tsd: Optional[Tsd] = None
    if best_tir is not None:
        flank_start      = max(1, best_tir.left_start - 50)
        flank_end_left   = best_tir.left_start - 1
        flank_start_right = best_tir.right_end + 1
        flank_end        = min(clen, best_tir.right_end + 50)
        left_flank  = (fa.fetch(contig, (flank_start, flank_end_left))
                       if flank_end_left >= flank_start else "")
        right_flank = (fa.fetch(contig, (flank_start_right, flank_end))
                       if flank_end >= flank_start_right else "")
        tsd = find_tsd(left_flank, right_flank,
                       cfg["tsd_min"], cfg["tsd_max"], cfg["tsd_max_slide"])

    # ── Stage 5: GC of PLV sequence ──────────────────────────────────────────
    plv_seq = fa.fetch(contig, (plv_start, plv_end))
    gc_plv  = gc_mean_of_seq(plv_seq)

    plv = Plv(
        plv_id="TBD",
        contig=contig, contig_length=clen,
        start=plv_start, end=plv_end, length=plv_length,
        tir=best_tir, tsd=tsd, orfs=plv_orfs,
        n_families=len(fams_in_plv),
        families_present=sorted(fams_in_plv),
        mcp_best_bitscore=mcp_best,
        gc_plv=gc_plv,
        has_tir=(best_tir is not None),
    )
    return dict(
        status="ok", cluster_index=ci, contig=contig,
        log_msgs=log_msgs, plv=plv,
    )


# ── Argument parsing ──────────────────────────────────────────────────────────
def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="findPLV.py",
        description="Identify Polinton-Like Viruses in eukaryotic genome assemblies.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Positional
    p.add_argument("genome", type=Path,
                   help="Input genome assembly FASTA (gzip OK)")

    # Database group
    db = p.add_argument_group("Database")
    db.add_argument("--mcp-plv", type=Path, required=True,
                    help="PLV major capsid protein HMM (required)")
    db.add_argument("--mcp-ncldv-virophage", type=Path, required=True,
                    help="NCLDV/Virophage MCP HMM for exclusion filtering (required)")
    db.add_argument("--pfam", type=Path, required=True,
                    help="Pfam-A.hmm database (required)")
    db.add_argument("--mcp-minor", type=Path, default=None,
                    help="PLV minor capsid protein HMM (optional)")

    # Optional
    opt = p.add_argument_group("Optional")
    opt.add_argument("--prefix", type=str, default="findPLV",
                     help="Output prefix for PLV IDs and FASTA headers  [default: findPLV]")
    opt.add_argument("--outdir", type=Path, default=Path("findPLV"),
                     help="Output directory                              [default: ./findPLV]")
    opt.add_argument("--seed-window", type=int, default=DEFAULTS["window_size"],
                     help="Seeding window size in bp centred on anchors  [default: 20000]")
    opt.add_argument("-t", "--threads", type=int, default=8,
                     help="CPU threads for HMM scanning and ORF prediction [default: 8]")
    opt.add_argument("--einverted-jobs", type=int, default=None,
                     help="Parallel einverted workers                    [default: --threads]")

    # Advanced — hidden from --help; adjust via DEFAULTS block at top of script
    p.add_argument("--min-contig",                 type=int,   default=DEFAULTS["min_contig"],                    help=argparse.SUPPRESS)
    p.add_argument("--min-plv-length",             type=int,   default=DEFAULTS["min_plv_length"],                help=argparse.SUPPRESS)
    p.add_argument("--min-families",               type=int,   default=DEFAULTS["min_families"],                  help=argparse.SUPPRESS)
    p.add_argument("--cluster-merge-gap",          type=int,   default=DEFAULTS["cluster_merge_gap"],             help=argparse.SUPPRESS)
    p.add_argument("--max-cluster-span",           type=int,   default=DEFAULTS["max_cluster_span"],              help=argparse.SUPPRESS)
    p.add_argument("--max-plv-length",             type=int,   default=DEFAULTS["max_plv_length"],                help=argparse.SUPPRESS)
    p.add_argument("--edge-gap-trim",              type=int,   default=DEFAULTS["edge_gap_trim"],                 help=argparse.SUPPRESS)
    p.add_argument("--gc-flank-search",            type=int,   default=DEFAULTS["gc_flank_search"],               help=argparse.SUPPRESS)
    p.add_argument("--gc-min-delta",               type=float, default=DEFAULTS["gc_min_delta_pct"],              help=argparse.SUPPRESS)
    p.add_argument("--mcp-evalue",                 type=float, default=DEFAULTS["mcp_evalue"],                    help=argparse.SUPPRESS)
    p.add_argument("--mcp-minor-evalue",           type=float, default=DEFAULTS["mcp_minor_evalue"],              help=argparse.SUPPRESS)
    p.add_argument("--pfam-evalue",                type=float, default=DEFAULTS["pfam_evalue"],                   help=argparse.SUPPRESS)
    p.add_argument("--ncldv-virophage-mcp-evalue", type=float, default=DEFAULTS["ncldv_virophage_mcp_evalue"],   help=argparse.SUPPRESS)
    p.add_argument("--dedup-overlap",              type=float, default=DEFAULTS["dedup_min_reciprocal_overlap"],  help=argparse.SUPPRESS)
    p.add_argument("--step",                       type=int,   default=DEFAULTS["window_step"],                   help=argparse.SUPPRESS)

    return p.parse_args(argv)


# ── Entry point ───────────────────────────────────────────────────────────────
def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    # Build runtime config from defaults + CLI overrides.
    cfg = dict(DEFAULTS)
    cfg["min_contig"]                 = args.min_contig
    cfg["min_plv_length"]             = args.min_plv_length
    cfg["window_size"]                = args.seed_window
    cfg["window_step"]                = args.step
    cfg["min_families"]               = args.min_families
    cfg["cluster_merge_gap"]          = args.cluster_merge_gap
    cfg["max_cluster_span"]           = args.max_cluster_span
    cfg["max_plv_length"]             = args.max_plv_length
    cfg["edge_gap_trim"]              = args.edge_gap_trim
    cfg["gc_flank_search"]            = args.gc_flank_search
    cfg["gc_min_delta_pct"]           = args.gc_min_delta
    cfg["mcp_evalue"]                 = args.mcp_evalue
    cfg["mcp_minor_evalue"]           = args.mcp_minor_evalue
    cfg["pfam_evalue"]                = args.pfam_evalue
    cfg["ncldv_virophage_mcp_evalue"] = args.ncldv_virophage_mcp_evalue

    if args.threads < 1:
        print("--threads must be >= 1", file=sys.stderr)
        return 2
    einverted_jobs = args.einverted_jobs if args.einverted_jobs is not None else args.threads
    if einverted_jobs < 1:
        print("--einverted-jobs must be >= 1", file=sys.stderr)
        return 2

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "tracks").mkdir(exist_ok=True)
    (args.outdir / "gggenomes").mkdir(exist_ok=True)

    setup_logging(args.outdir / "run.log")
    t0 = time.time()
    log_info(
        f"findPLV started | prefix='{args.prefix}' | "
        f"threads={args.threads} | einverted_jobs={einverted_jobs} | "
        f"input={args.genome}"
    )
    log_info(
        f"Parameters | min_contig={cfg['min_contig']:,} bp | "
        f"min_plv_length={cfg['min_plv_length']:,} bp | "
        f"window={cfg['window_size']:,} bp | "
        f"min_families={cfg['min_families']} | "
        f"cluster_merge_gap={cfg['cluster_merge_gap']:,} bp | "
        f"max_cluster_span={cfg['max_cluster_span']:,} bp | "
        f"max_plv_length={cfg['max_plv_length']:,} bp | "
        f"edge_gap_trim={cfg['edge_gap_trim']:,} bp | "
        f"gc_flank_search={cfg['gc_flank_search']:,} bp | "
        f"gc_min_delta={cfg['gc_min_delta_pct']:.1f} pp | "
        f"mcp_evalue={cfg['mcp_evalue']:g} | "
        f"pfam_evalue={cfg['pfam_evalue']:g}"
    )

    # File existence checks
    for path, kind in [
        (args.genome,               "genome"),
        (args.mcp_plv,              "MCP HMM"),
        (args.mcp_ncldv_virophage,  "NCLDV/Virophage MCP HMM"),
        (args.pfam,                 "Pfam-A HMM"),
    ]:
        if not path.exists():
            log_error(f"Missing {kind} file: {path}")
            return 2
    if args.mcp_minor is not None and not args.mcp_minor.exists():
        log_error(f"Missing mCP HMM file: {args.mcp_minor}")
        return 2
    if not shutil.which("einverted"):
        log_error(
            "einverted not found in PATH. "
            "Install EMBOSS (conda install -c bioconda emboss)."
        )
        return 2

    family_table = build_family_table()

    # Anchor families: base set + mCP if the HMM was supplied.
    anchor_families = set(ANCHOR_FAMILIES_BASE)
    if args.mcp_minor is not None:
        anchor_families.add(MCP_MINOR_LABEL)
        log_info(
            f"mCP HMM supplied ({args.mcp_minor.name}); "
            f"mCP added to anchor families and seeding"
        )
    else:
        log_info(
            "No --mcp-minor HMM supplied; mCP scanning disabled. "
            "Anchor families: " + ", ".join(sorted(anchor_families))
        )
    anchor_families = frozenset(anchor_families)
    cfg["anchor_families"] = anchor_families

    # ── Stage 1: ORF prediction ───────────────────────────────────────────────
    orfs_by_id, contig_lengths = predict_orfs(args.genome, cfg["min_contig"], args.threads)

    # ── Stage 2a: major capsid scan ──────────────────────────────────────────
    scan_mcp(orfs_by_id, args.mcp_plv, cfg["mcp_evalue"], args.threads)

    # ── Stage 2b: minor capsid scan (optional) ────────────────────────────────
    if args.mcp_minor is not None:
        scan_mcp_minor(orfs_by_id, args.mcp_minor, cfg["mcp_minor_evalue"], args.threads)

    # Select contigs with at least one capsid hit for the Pfam scan.
    seed_contigs = {
        o.contig for o in orfs_by_id.values()
        if o.family in (MCP_LABEL, MCP_MINOR_LABEL)
    }
    if not seed_contigs:
        log_warning("No MCP or mCP hits detected. No PLV candidates to report.")
        _write_empty_outputs(args, t0)
        return 0
    log_info(f"Capsid-positive contigs selected for Pfam scan: {len(seed_contigs):,}")

    # ── Stage 2c: scoped Pfam scan on capsid-positive contigs ────────────────
    pfam_targets = [o for o in orfs_by_id.values() if o.contig in seed_contigs]
    scan_pfam(pfam_targets, args.pfam, family_table, cfg["pfam_evalue"], args.threads)

    # ── Stage 2d: multi-anchor seeding ───────────────────────────────────────
    clusters = find_seed_regions(
        orfs_by_id         = orfs_by_id,
        anchor_families    = anchor_families,
        window_size        = cfg["window_size"],
        min_families       = cfg["min_families"],
        require_mcp        = cfg["require_mcp"],
        cluster_merge_gap  = cfg["cluster_merge_gap"],
        max_cluster_span   = cfg["max_cluster_span"],
    )
    if not clusters:
        log_warning("No marker clusters passed seeding criteria. No PLV candidates to report.")
        _write_empty_outputs(args, t0)
        return 0

    # ── Stage 2e: NCLDV / Virophage MCP exclusion ────────────────────────────
    clusters = filter_clusters_by_ncldv_virophage(
        clusters, orfs_by_id,
        args.mcp_ncldv_virophage,
        cfg["ncldv_virophage_mcp_evalue"],
        args.threads,
    )
    if not clusters:
        log_warning(
            "All seed clusters excluded by NCLDV/Virophage MCP filter. "
            "No PLV candidates to report."
        )
        _write_empty_outputs(args, t0)
        return 0

    # ── Stages 3-5: TIR / TSD / GC per cluster (parallel) ───────────────────
    fa = pyfastx.Fasta(str(args.genome), build_index=True, uppercase=True)
    orfs_by_contig: Dict[str, List[Orf]] = {}
    for o in orfs_by_id.values():
        orfs_by_contig.setdefault(o.contig, []).append(o)

    tasks: List[dict] = []
    for ci, cl in enumerate(clusters, start=1):
        tasks.append(dict(
            cluster_index = ci,
            contig        = cl["contig"],
            cluster_start = cl["cluster_start"],
            cluster_end   = cl["cluster_end"],
            contig_length = contig_lengths[cl["contig"]],
            cfg           = cfg,
            genome_path   = str(args.genome),
            contig_orfs   = orfs_by_contig.get(cl["contig"], []),
        ))

    plvs: List[Plv] = []
    n_workers = max(1, min(einverted_jobs, len(tasks)))
    log_info(
        f"TIR detection: dispatching {len(tasks)} cluster(s) "
        f"across {n_workers} parallel einverted worker(s)"
    )

    if n_workers == 1:
        results_iter = (_process_cluster(t) for t in tasks)
    else:
        executor = ProcessPoolExecutor(max_workers=n_workers)
        results_iter = executor.map(_process_cluster, tasks, chunksize=1)

    n_done = 0
    fatal_msg: Optional[str] = None
    try:
        for res in results_iter:
            n_done += 1
            for level, msg in res.get("log_msgs", []) or []:
                (log_warning if level == "warning" else log_info)(msg)
            status = res.get("status")
            if status == "fatal":
                fatal_msg = res.get("message", "unknown fatal error in worker")
                log_error(fatal_msg)
                break
            if status == "skip":
                log_warning(res.get("message", f"Cluster {n_done}: discarded"))
                continue
            if status == "ok":
                plv = res["plv"]
                plvs.append(plv)
                tir_str = (
                    f"TIR={plv.tir.tir_length} bp @ {plv.tir.tir_identity:.1f}% id"
                    if plv.tir else "TIR=none"
                )
                log_info(
                    f"Cluster {res['cluster_index']}/{len(tasks)} accepted | "
                    f"{plv.contig}:{plv.start:,}-{plv.end:,} | "
                    f"length={plv.length:,} bp | {tir_str} | "
                    f"marker_families={plv.n_families} | "
                    f"TSD={'yes' if plv.tsd else 'no'}"
                )
    finally:
        if n_workers > 1:
            executor.shutdown(wait=True)

    if fatal_msg is not None:
        return 2

    # ── Deduplication ─────────────────────────────────────────────────────────
    plvs = deduplicate_plvs(plvs, args.dedup_overlap)

    # ── Sort: TIR-bearing first, then by genome position ─────────────────────
    plvs.sort(key=lambda p: (not p.has_tir, p.contig, p.start, p.end))
    for i, p in enumerate(plvs, start=1):
        p.plv_id = f"{args.prefix}_PLV_{i:03d}"

    # Warn on remaining overlaps (rare; post-dedup).
    plvs_by_contig: Dict[str, List[Plv]] = {}
    for p in plvs:
        plvs_by_contig.setdefault(p.contig, []).append(p)
    for contig, lst in plvs_by_contig.items():
        lst.sort(key=lambda x: x.start)
        for a, b in zip(lst, lst[1:]):
            if b.start <= a.end:
                log_warning(
                    f"Overlapping PLV spans on {contig}: "
                    f"{a.plv_id} ({a.start:,}-{a.end:,}) overlaps "
                    f"{b.plv_id} ({b.start:,}-{b.end:,})"
                )

    n_with_tir = sum(1 for p in plvs if p.has_tir)
    log_info(
        f"Final: {len(plvs)} PLV candidate(s) from {len(tasks)} cluster(s) "
        f"({n_with_tir} with TIR, {len(plvs) - n_with_tir} without TIR)"
    )

    # ── Stage 6: write outputs ────────────────────────────────────────────────
    results_path  = args.outdir / "plv_results.tsv"
    summary_path  = args.outdir / "run_summary.txt"
    fasta_path    = args.outdir / "plvs.fna"
    gff3_path     = args.outdir / "plvs.gff3"
    seqs_gg_path  = args.outdir / "gggenomes" / "seqs.tsv"
    genes_gg_path = args.outdir / "gggenomes" / "genes.tsv"
    bedgraph_path = args.outdir / "tracks" / "gc_500bp.bedgraph"
    log_path      = args.outdir / "run.log"

    write_plv_tsv(plvs, fa, results_path, args.prefix)
    write_multifasta(plvs, fa, fasta_path, args.prefix)
    write_gff3(plvs, gff3_path)
    write_gggenomes(plvs, seqs_gg_path, genes_gg_path)
    write_bedgraph(plvs, fa, contig_lengths, bedgraph_path,
                   gc_window=cfg["gc_window"], flank=cfg["gc_flank"])

    elapsed = time.time() - t0
    write_run_summary_txt(plvs, summary_path, args.prefix, args.genome, elapsed)

    log_output(f"Results table    -> {results_path}")
    log_output(f"Run summary      -> {summary_path}")
    log_output(f"PLV sequences    -> {fasta_path}")
    log_output(f"Annotations GFF3 -> {gff3_path}")
    log_output(f"gggenomes seqs   -> {seqs_gg_path}")
    log_output(f"gggenomes genes  -> {genes_gg_path}")
    log_output(f"GC bedgraph      -> {bedgraph_path}")
    log_output(f"Run log          -> {log_path}")

    timing = (
        f"{elapsed / 3600:.2f} h ({elapsed:.1f} s)" if elapsed >= 3600
        else f"{elapsed / 60:.2f} min ({elapsed:.1f} s)" if elapsed >= 60
        else f"{elapsed:.1f} s"
    )
    log_info(f"findPLV completed in {timing}")
    return 0


def _write_empty_outputs(args: argparse.Namespace, t0: float) -> None:
    """Write empty output files and a summary when no PLV candidates are found."""
    write_plv_tsv([], None, args.outdir / "plv_results.tsv", args.prefix)
    write_run_summary_txt([], args.outdir / "run_summary.txt",
                          args.prefix, args.genome, time.time() - t0)


if __name__ == "__main__":
    sys.exit(main())
