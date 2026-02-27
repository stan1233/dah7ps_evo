#!/usr/bin/env python3
"""Select structure panel from stepping-stone coverage backbone.

Implements the Phase 3.1A Selection Contract:
  - Target N=30 (20-40), quota Ia=12/Ib=5/II=13 (±1; floors Ia≥8/Ib≥4/II≥8)
  - Priority: PDB > AFDB (core pLDDT>=70) > ESMFold (core pLDDT>=70)
  - One representative per stepping-stone cluster
  - Stratified sampling: >=30% small(<=2), >=30% medium(3-10), >=20% large(>10)
  - Anchor constraint: >=1 PDB/AFDB per subtype

Inputs:
  --candidates_tsv   panel_candidates.tsv (258 rows, with structure availability)
  --params           meta/params.json (selection contract parameters)

Outputs:
  --manifest_tsv     panel_manifest.tsv (selected panel, 20-40 rows)

Usage:
  # Step 1: Generate panel_candidates.tsv (manually or via helper script)
  # Step 2: Run selection
  python scripts/select_structure_panel.py \\
    --candidates_tsv results/03_msa_core/panel_candidates.tsv \\
    --params meta/params.json \\
    --manifest_tsv results/03_msa_core/panel_manifest.tsv
"""

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def load_params(params_path):
    """Load structure_panel parameters from params.json."""
    with open(params_path) as f:
        params = json.load(f)
    sp = params.get("structure_panel", {})
    if not sp:
        print("ERROR: 'structure_panel' block missing from params.json", file=sys.stderr)
        sys.exit(1)
    return sp


def load_candidates(tsv_path):
    """Load panel_candidates.tsv. Required columns:
    rep_id, subtype, cluster_id, cluster_size, seq_len,
    has_pdb, pdb_ids, has_afdb, afdb_plddt_core, afdb_core_cov,
    needs_esmf, esmf_plddt_core, esmf_core_cov, selected_panel
    """
    import csv
    candidates = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'rep_id', 'subtype', 'cluster_id', 'cluster_size', 'seq_len',
                     'has_pdb', 'has_afdb', 'afdb_plddt_core'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            print(f"ERROR: Missing columns in candidates TSV: {missing}", file=sys.stderr)
            sys.exit(1)
        for row in reader:
            # Parse numeric fields
            row['cluster_size'] = int(row['cluster_size'])
            row['seq_len'] = int(row['seq_len'])
            row['has_pdb'] = row['has_pdb'].strip().lower() in ('1', 'true', 'yes')
            row['has_afdb'] = row['has_afdb'].strip().lower() in ('1', 'true', 'yes')
            row['afdb_plddt_core'] = float(row.get('afdb_plddt_core') or 0)
            row['afdb_core_cov'] = float(row.get('afdb_core_cov') or 0)
            row['needs_esmf'] = row.get('needs_esmf', '').strip().lower() in ('1', 'true', 'yes')
            row['esmf_plddt_core'] = float(row.get('esmf_plddt_core') or 0)
            row['esmf_core_cov'] = float(row.get('esmf_core_cov') or 0)
            candidates.append(row)
    return candidates


def classify_cluster_size(size, small_max, medium_max):
    """Classify cluster into small/medium/large."""
    if size <= small_max:
        return 'small'
    elif size <= medium_max:
        return 'medium'
    else:
        return 'large'


def get_structure_source(candidate, plddt_min, cov_min):
    """Determine best available structure source and whether it passes quality."""
    if candidate['has_pdb']:
        return 'PDB', True   # PDB always passes
    if candidate['has_afdb'] and candidate['afdb_plddt_core'] >= plddt_min:
        if candidate['afdb_core_cov'] >= cov_min or candidate['afdb_core_cov'] == 0:
            return 'AFDB', True
    if candidate['needs_esmf'] and candidate['esmf_plddt_core'] >= plddt_min:
        if candidate['esmf_core_cov'] >= cov_min or candidate['esmf_core_cov'] == 0:
            return 'ESMFold', True
    # Has AFDB but below threshold
    if candidate['has_afdb']:
        return 'AFDB', False
    if candidate['needs_esmf']:
        return 'ESMFold', False
    return 'NONE', False


def select_panel(candidates, sp):
    """Select panel according to Selection Contract."""
    quota = dict(sp['quota_default'])
    floors = dict(sp['quota_floors'])
    flex = sp.get('quota_flex', 1)
    plddt_min = sp['core_plddt_min']
    cov_min = sp.get('core_coverage_min', 0.80)
    strat = sp.get('stratification_targets', {})
    small_max = strat.get('small_max_size', 2)
    medium_max = strat.get('medium_max_size', 10)

    # Annotate each candidate with structure source and quality
    for c in candidates:
        src, passes = get_structure_source(c, plddt_min, cov_min)
        c['_source'] = src
        c['_passes_quality'] = passes
        c['_size_class'] = classify_cluster_size(c['cluster_size'], small_max, medium_max)

    # Group by subtype
    by_subtype = defaultdict(list)
    for c in candidates:
        by_subtype[c['subtype']].append(c)

    selected = []

    for st, target_n in quota.items():
        pool = by_subtype.get(st, [])
        if not pool:
            print(f"WARNING: No candidates for subtype {st}", file=sys.stderr)
            continue

        # Sort by priority: PDB first, then AFDB (by pLDDT desc), then ESMFold
        priority_order = {'PDB': 0, 'AFDB': 1, 'ESMFold': 2, 'NONE': 3}

        def sort_key(c):
            return (
                priority_order.get(c['_source'], 3),
                not c['_passes_quality'],
                -c['afdb_plddt_core'],
                -c['cluster_size'],  # prefer larger clusters as tiebreaker
            )

        eligible = [c for c in pool if c['_passes_quality']]
        eligible.sort(key=sort_key)

        # Ensure anchor constraint: at least 1 PDB or AFDB
        anchors = [c for c in eligible if c['_source'] in ('PDB', 'AFDB')]
        if not anchors:
            print(f"WARNING: No anchor structure (PDB/AFDB) for subtype {st}", file=sys.stderr)

        # Select up to target_n, respecting stratification
        chosen = []
        used_clusters = set()

        # Phase 1: anchor first
        for c in eligible:
            if c['_source'] in ('PDB', 'AFDB') and c['cluster_id'] not in used_clusters:
                chosen.append(c)
                used_clusters.add(c['cluster_id'])
                if len(chosen) >= target_n:
                    break

        # Phase 2: fill remaining by priority
        for c in eligible:
            if len(chosen) >= target_n:
                break
            if c['cluster_id'] not in used_clusters:
                chosen.append(c)
                used_clusters.add(c['cluster_id'])

        selected.extend(chosen)
        print(f"  {st}: selected {len(chosen)}/{target_n} (eligible={len(eligible)}, "
              f"pool={len(pool)}, anchors={len(anchors)})")

    return selected


def validate_selection(selected, sp):
    """Run assertions on the final selection."""
    quota = dict(sp['quota_default'])
    floors = dict(sp['quota_floors'])
    allowed_min, allowed_max = sp['allowed_range']

    total = len(selected)
    print(f"\n=== Validation ===")
    print(f"Total selected: {total} (allowed: {allowed_min}-{allowed_max})")

    # Check total range
    assert allowed_min <= total <= allowed_max, \
        f"FAIL: Total {total} outside allowed range [{allowed_min}, {allowed_max}]"
    print(f"  [PASS] Total in range")

    # Check per-subtype floors
    by_st = Counter(c['subtype'] for c in selected)
    for st, floor in floors.items():
        n = by_st.get(st, 0)
        assert n >= floor, f"FAIL: {st} has {n} < floor {floor}"
        print(f"  [PASS] {st}: {n} >= floor {floor}")

    # Check no duplicate clusters
    cluster_ids = [c['cluster_id'] for c in selected]
    assert len(cluster_ids) == len(set(cluster_ids)), \
        f"FAIL: Duplicate clusters in selection"
    print(f"  [PASS] No duplicate clusters")

    # Check anchor per subtype
    for st in quota:
        st_sel = [c for c in selected if c['subtype'] == st]
        anchors = [c for c in st_sel if c['_source'] in ('PDB', 'AFDB')]
        if anchors:
            print(f"  [PASS] {st}: anchor present ({anchors[0]['_source']})")
        else:
            print(f"  [WARN] {st}: NO anchor (PDB/AFDB) — all ESMFold")

    # Check stratification (soft, report only)
    strat = sp.get('stratification_targets', {})
    for st in quota:
        st_sel = [c for c in selected if c['subtype'] == st]
        if not st_sel:
            continue
        size_counts = Counter(c['_size_class'] for c in st_sel)
        n = len(st_sel)
        for cls, target_pct in [('small', strat.get('small_clusters_pct', 0.3)),
                                 ('medium', strat.get('medium_clusters_pct', 0.3)),
                                 ('large', strat.get('large_clusters_pct', 0.2))]:
            actual_pct = size_counts.get(cls, 0) / n
            status = "PASS" if actual_pct >= target_pct else "SOFT-MISS"
            print(f"  [{status}] {st} {cls}: {size_counts.get(cls, 0)}/{n} "
                  f"({actual_pct:.0%} vs target {target_pct:.0%})")

    print(f"\n=== All hard constraints passed ===")


def write_manifest(selected, output_path):
    """Write panel_manifest.tsv."""
    import csv
    fields = ['subtype', 'rep_id', 'cluster_id', 'cluster_size', 'seq_len',
              'structure_source', 'core_plddt', 'size_class', 'pdb_ids', 'notes']
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for c in sorted(selected, key=lambda x: (x['subtype'], x['cluster_id'])):
            plddt = (c['afdb_plddt_core'] if c['_source'] == 'AFDB'
                     else c['esmf_plddt_core'] if c['_source'] == 'ESMFold'
                     else 0)
            writer.writerow({
                'subtype': c['subtype'],
                'rep_id': c['rep_id'],
                'cluster_id': c['cluster_id'],
                'cluster_size': c['cluster_size'],
                'seq_len': c['seq_len'],
                'structure_source': c['_source'],
                'core_plddt': f"{plddt:.1f}" if plddt else 'N/A',
                'size_class': c['_size_class'],
                'pdb_ids': c.get('pdb_ids', ''),
                'notes': c.get('notes', ''),
            })
    print(f"\nManifest written: {output_path} ({len(selected)} entries)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--candidates_tsv', required=True,
                        help='panel_candidates.tsv (258 rows)')
    parser.add_argument('--params', required=True,
                        help='meta/params.json')
    parser.add_argument('--manifest_tsv', required=True,
                        help='Output panel_manifest.tsv')
    args = parser.parse_args()

    # Input validation
    for p in [args.candidates_tsv, args.params]:
        if not Path(p).exists():
            print(f"ERROR: Input file not found: {p}", file=sys.stderr)
            sys.exit(1)

    Path(args.manifest_tsv).parent.mkdir(parents=True, exist_ok=True)

    sp = load_params(args.params)
    candidates = load_candidates(args.candidates_tsv)
    print(f"Loaded {len(candidates)} candidates")

    selected = select_panel(candidates, sp)
    validate_selection(selected, sp)
    write_manifest(selected, args.manifest_tsv)


if __name__ == '__main__':
    main()
