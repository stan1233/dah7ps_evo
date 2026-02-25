#!/usr/bin/env python3
"""
HMMER Utilities — Reusable functions for domtblout parsing and hit stitching.

Used by:
  - Phase 2.1: qc_length_coverage.py (merged HMM coverage for Type II)
  - Phase 3.6: extract_core_domains.py (hit stitching for core domain extraction)

Key concepts:
  - "Hit stitching" merges multiple domain hits from the same query sequence
    into a single coverage region, solving the α2β3 insertion fragmentation
    problem in Type II DAH7PS (SOP CHECK-06).
  - Coverage is computed in HMM coordinate space (hmm_from/hmm_to), not
    sequence coordinate space, to correctly measure how much of the model
    is covered.
"""


def parse_domtbl_all_domains(path, ievalue_max=1e-5, min_hmm_span=30):
    """Parse domtblout, return all qualifying domain hits per sequence.

    Args:
        path: Path to HMMER --domtblout output file.
        ievalue_max: Maximum individual domain E-value to retain (default 1e-5).
        min_hmm_span: Minimum HMM span (hmm_to - hmm_from + 1) to retain (default 30).

    Returns:
        dict: seqid -> list of domain dicts, each with keys:
              full_score, dom_score, i_evalue,
              hmm_from, hmm_to, env_from, env_to
    """
    domains = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split(maxsplit=22)
            if len(parts) < 22:
                continue

            seqid = parts[0]
            full_score = float(parts[7])
            i_evalue = float(parts[12])
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            env_from = int(parts[19])
            env_to = int(parts[20])

            # Filter: skip weak/tiny hits
            hmm_span = hmm_to - hmm_from + 1
            if i_evalue > ievalue_max:
                continue
            if hmm_span < min_hmm_span:
                continue

            dom = {
                "full_score": full_score,
                "dom_score": dom_score,
                "i_evalue": i_evalue,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "env_from": env_from,
                "env_to": env_to,
            }

            if seqid not in domains:
                domains[seqid] = []
            domains[seqid].append(dom)

    return domains


def parse_domtbl_best_domain(path):
    """Parse domtblout, return best single domain hit per sequence.

    Returns:
        dict: seqid -> (dom_score, hmm_from, hmm_to, env_from, env_to)
    """
    best = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split(maxsplit=22)
            if len(parts) < 22:
                continue

            seqid = parts[0]
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            env_from = int(parts[19])
            env_to = int(parts[20])

            if seqid not in best or dom_score > best[seqid][0]:
                best[seqid] = (dom_score, hmm_from, hmm_to, env_from, env_to)
    return best


def merge_intervals(intervals, merge_gap=0):
    """Merge overlapping/adjacent intervals.

    Args:
        intervals: List of (start, end) tuples.
        merge_gap: Maximum gap between intervals to still merge (default 0).

    Returns:
        List of [start, end] merged intervals, sorted.
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = []
    for s, e in intervals:
        if not merged or s > merged[-1][1] + merge_gap + 1:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return merged


def union_coverage_hmm(domains, hmm_len, merge_gap=0):
    """Compute union HMM coverage from multiple domain hits.

    Args:
        domains: List of domain dicts (must have 'hmm_from', 'hmm_to').
        hmm_len: Length of the HMM model (number of match states).
        merge_gap: Gap tolerance for merging adjacent intervals (default 0).

    Returns:
        (coverage_fraction, merged_intervals, n_domains_used)
    """
    if not domains:
        return 0.0, [], 0

    intervals = [(d["hmm_from"], d["hmm_to"]) for d in domains]
    merged = merge_intervals(intervals, merge_gap)
    cov_len = sum(e - s + 1 for s, e in merged)
    cov_len = min(cov_len, hmm_len)  # cap at model length

    return cov_len / hmm_len, merged, len(domains)


def union_coverage_env(domains, merge_gap=0):
    """Compute union envelope coverage on the sequence coordinate.

    Useful for Phase 3.6 core domain extraction (stitched env_from..env_to).

    Args:
        domains: List of domain dicts (must have 'env_from', 'env_to').
        merge_gap: Gap tolerance for merging (default 0).

    Returns:
        (total_residues, merged_intervals)
    """
    if not domains:
        return 0, []

    intervals = [(d["env_from"], d["env_to"]) for d in domains]
    merged = merge_intervals(intervals, merge_gap)
    total = sum(e - s + 1 for s, e in merged)

    return total, merged
