#!/usr/bin/env python3
"""Build panel_candidates.tsv for Phase 3.1A-1.

Steps:
  1. Build skeleton from local data (rep_id, subtype, cluster_id, cluster_size, seq_len)
  2. Query PDB availability via PDBe SIFTS API (UniProt -> PDB mapping)
  3. Query AFDB availability via AlphaFold DB API
  4. Compute core-region pLDDT from AFDB confidence (B-factor in mmCIF)
  5. Output panel_candidates.tsv

Usage:
  python3 scripts/build_panel_candidates.py --workdir /home/tynan/0218
"""

import argparse
import csv
import json
import os
import sys
import time
import urllib.request
import urllib.error
from collections import defaultdict
from pathlib import Path


def load_ids_from_fasta(fasta_path):
    """Extract (id, seq_len) from FASTA."""
    results = {}
    current_id = None
    current_len = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    results[current_id] = current_len
                current_id = line[1:].strip().split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_id:
        results[current_id] = current_len
    return results


def build_subtype_map(seeds60_files):
    """Build rep_id -> subtype from seeds60 FASTA files."""
    subtype_map = {}
    for st, path in seeds60_files:
        with open(path) as f:
            for line in f:
                if line.startswith('>'):
                    sid = line[1:].strip().split()[0]
                    subtype_map[sid] = st
    return subtype_map


def build_cluster_sizes(cluster_tsv):
    """Build rep_id -> cluster_size from stepping_stones_cluster.tsv."""
    sizes = defaultdict(int)
    with open(cluster_tsv) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                rep, _ = parts
                sizes[rep] += 1
    return dict(sizes)


def extract_accession(rep_id):
    """Extract UniProt accession from UniRef90 ID.
    UniRef90_A0A9Q5I2C3 -> A0A9Q5I2C3
    UniRef90_UPI003CF8EA9B -> UPI003CF8EA9B (UniParc, no AFDB)
    """
    acc = rep_id
    if acc.startswith('UniRef90_'):
        acc = acc[len('UniRef90_'):]
    return acc


def is_uniprot_accession(acc):
    """Check if accession looks like a UniProt ID (not UniParc UPI)."""
    return not acc.startswith('UPI')


def query_pdb_sifts(accession, cache_dir):
    """Query PDBe SIFTS for UniProt -> PDB mapping. Returns list of PDB IDs."""
    cache_file = cache_dir / f"{accession}_pdb.json"
    if cache_file.exists():
        with open(cache_file) as f:
            data = json.load(f)
        return data.get('pdb_ids', [])

    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{accession}"
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        # Extract PDB IDs
        pdb_ids = []
        for uniprot_acc, mapping in data.items():
            for pdb_entry in mapping.get('PDB', {}):
                pdb_ids.append(pdb_entry.upper())
        result = {'pdb_ids': list(set(pdb_ids)), 'status': 'found'}
    except urllib.error.HTTPError as e:
        if e.code == 404:
            result = {'pdb_ids': [], 'status': 'not_found'}
        else:
            result = {'pdb_ids': [], 'status': f'error_{e.code}'}
    except Exception as e:
        result = {'pdb_ids': [], 'status': f'error_{str(e)[:50]}'}

    with open(cache_file, 'w') as f:
        json.dump(result, f)
    return result.get('pdb_ids', [])


def query_afdb(accession, cache_dir):
    """Query AlphaFold DB for structure availability. Returns entry info or None."""
    cache_file = cache_dir / f"{accession}_afdb.json"
    if cache_file.exists():
        with open(cache_file) as f:
            data = json.load(f)
        return data if data.get('status') != 'not_found' else None

    url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        # data is a list of entries; take the first (latest version)
        if isinstance(data, list) and data:
            entry = data[0]
            result = {
                'status': 'found',
                'entryId': entry.get('entryId', ''),
                'cifUrl': entry.get('cifUrl', ''),
                'pdbUrl': entry.get('pdbUrl', ''),
                'paeImageUrl': entry.get('paeImageUrl', ''),
                'uniprotAccession': entry.get('uniprotAccession', accession),
                'uniprotId': entry.get('uniprotId', ''),
                'gene': entry.get('gene', ''),
                'organism': entry.get('organismScientificName', ''),
                'globalPlddt': entry.get('globalMetricValue', 0),
            }
        else:
            result = {'status': 'not_found'}
    except urllib.error.HTTPError as e:
        if e.code == 404 or e.code == 422:
            result = {'status': 'not_found'}
        else:
            result = {'status': f'error_{e.code}'}
    except Exception as e:
        result = {'status': f'error_{str(e)[:50]}'}

    with open(cache_file, 'w') as f:
        json.dump(result, f)
    return result if result.get('status') == 'found' else None


def get_hmm_core_coords(rep_id, domtbl_files, hmm_lengths):
    """Get HMM core region coordinates (sequence coords) from domtblout.
    Returns (env_from, env_to) for the best domain hit, or merged for Type II.
    """
    acc_raw = rep_id
    best_coords = None
    best_score = -1

    for st, domtbl_path, hmm_len in domtbl_files:
        if not os.path.exists(domtbl_path):
            continue
        with open(domtbl_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 22:
                    continue
                target_name = parts[0]
                if target_name != acc_raw:
                    continue
                score = float(parts[13])  # domain score
                ali_from = int(parts[17])  # ali coord start (1-based)
                ali_to = int(parts[18])    # ali coord end
                if score > best_score:
                    best_score = score
                    best_coords = (ali_from, ali_to)

    return best_coords


def download_afdb_pdb(accession, afdb_entry, struct_cache_dir):
    """Download AFDB PDB file for pLDDT extraction. Returns local path or None."""
    pdb_url = afdb_entry.get('pdbUrl', '')
    if not pdb_url:
        return None

    local_path = struct_cache_dir / f"AF-{accession}-F1-model_v4.pdb"
    if local_path.exists():
        return local_path

    try:
        urllib.request.urlretrieve(pdb_url, str(local_path))
        return local_path
    except Exception as e:
        print(f"  WARNING: Failed to download AFDB PDB for {accession}: {e}", file=sys.stderr)
        return None


def compute_core_plddt(pdb_path, core_from, core_to):
    """Compute mean pLDDT for core region from AFDB PDB file.
    AFDB stores pLDDT in the B-factor column.
    core_from/core_to are 1-based sequence coordinates.
    """
    if not pdb_path or not os.path.exists(pdb_path):
        return 0.0, 0.0

    residue_plddt = {}  # residue_num -> pLDDT (from CA atom)
    total_residues = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                resnum = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
                residue_plddt[resnum] = bfactor
                total_residues = max(total_residues, resnum)

    if not residue_plddt:
        return 0.0, 0.0

    # Extract core region pLDDT
    core_plddts = []
    for i in range(core_from, core_to + 1):
        if i in residue_plddt:
            core_plddts.append(residue_plddt[i])

    if not core_plddts:
        return 0.0, 0.0

    core_mean = sum(core_plddts) / len(core_plddts)
    core_coverage = len(core_plddts) / (core_to - core_from + 1)

    return core_mean, core_coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--workdir', default='/home/tynan/0218')
    parser.add_argument('--rate_limit', type=float, default=0.2,
                        help='Seconds between API calls (default: 0.2)')
    parser.add_argument('--skip_download', action='store_true',
                        help='Skip AFDB PDB downloads (use cached or leave pLDDT=0)')
    args = parser.parse_args()

    wd = Path(args.workdir)
    outdir = wd / 'results' / '03_msa_core'
    outdir.mkdir(parents=True, exist_ok=True)

    # Cache dirs for API responses and structure files
    cache_dir = outdir / 'structure_availability' / 'raw'
    cache_dir.mkdir(parents=True, exist_ok=True)
    struct_cache_dir = outdir / 'structure_availability' / 'afdb_pdb'
    struct_cache_dir.mkdir(parents=True, exist_ok=True)

    # === Step A: Build skeleton table ===
    print("=== Step A: Building skeleton table ===")

    # Load rep sequences
    reps = load_ids_from_fasta(wd / 'results/02_qc/stepping_stones_rep_seq.fasta')
    print(f"  Loaded {len(reps)} reps")

    # Load subtype map
    subtype_map = build_subtype_map([
        ('Ia', wd / 'results/02_qc/seeds60_Ia.fasta'),
        ('Ib', wd / 'results/02_qc/seeds60_Ib.fasta'),
        ('II', wd / 'results/02_qc/seeds60_II.fasta'),
    ])

    # Load cluster sizes
    cluster_sizes = build_cluster_sizes(wd / 'results/02_qc/stepping_stones_cluster.tsv')

    # Build skeleton rows
    rows = []
    for rep_id, seq_len in reps.items():
        rows.append({
            'rep_id': rep_id,
            'subtype': subtype_map.get(rep_id, 'Unknown'),
            'cluster_id': rep_id,  # rep is its own cluster ID
            'cluster_size': cluster_sizes.get(rep_id, 1),
            'seq_len': seq_len,
        })

    print(f"  Skeleton: {len(rows)} rows, subtypes: "
          + ", ".join(f"{st}={sum(1 for r in rows if r['subtype']==st)}" for st in ['Ia','Ib','II']))

    # === Step B: Query PDB and AFDB availability ===
    print("\n=== Step B: Querying structure availability ===")

    # Determine which domtblout files to use for core coords
    domtbl_files = [
        ('Ia', str(wd / 'results/01_mining/hits_Ia.domtbl'), 355),
        ('Ib', str(wd / 'results/01_mining/hits_Ib_vs_dah7ps_v41.domtbl'), 334),
        ('II', str(wd / 'results/01_mining/hits_II.domtbl'), 471),
    ]
    hmm_lengths = {'Ia': 355, 'Ib': 334, 'II': 471}

    total = len(rows)
    for i, row in enumerate(rows):
        acc = extract_accession(row['rep_id'])
        is_uniprot = is_uniprot_accession(acc)

        # Progress
        if (i + 1) % 20 == 0 or i == 0:
            print(f"  [{i+1}/{total}] Processing {acc} ({row['subtype']})...")

        # PDB query
        pdb_ids = query_pdb_sifts(acc, cache_dir) if is_uniprot else []
        row['has_pdb'] = 1 if pdb_ids else 0
        row['pdb_ids'] = ','.join(pdb_ids[:5]) if pdb_ids else ''
        if pdb_ids:
            time.sleep(args.rate_limit)

        # AFDB query
        if is_uniprot:
            afdb_entry = query_afdb(acc, cache_dir)
        else:
            afdb_entry = None

        row['has_afdb'] = 1 if afdb_entry else 0
        row['afdb_global_plddt'] = afdb_entry.get('globalPlddt', 0) if afdb_entry else 0

        # Core region coordinates (from domtblout)
        st = row['subtype']
        st_domtbls = [(s, d, h) for s, d, h in domtbl_files if s == st]
        core_coords = get_hmm_core_coords(row['rep_id'], st_domtbls, hmm_lengths)

        # Download AFDB PDB and compute core pLDDT
        if afdb_entry and not args.skip_download and core_coords:
            pdb_path = download_afdb_pdb(acc, afdb_entry, struct_cache_dir)
            core_from, core_to = core_coords
            core_plddt, core_cov = compute_core_plddt(pdb_path, core_from, core_to)
            row['afdb_plddt_core'] = round(core_plddt, 1)
            row['afdb_core_cov'] = round(core_cov, 3)
        elif afdb_entry and args.skip_download:
            row['afdb_plddt_core'] = row['afdb_global_plddt']  # fallback to global
            row['afdb_core_cov'] = 0  # unknown
        else:
            row['afdb_plddt_core'] = 0
            row['afdb_core_cov'] = 0

        # Determine needs_esmf
        if row['has_pdb']:
            row['needs_esmf'] = 0
        elif row['has_afdb'] and row['afdb_plddt_core'] >= 70 and (row['afdb_core_cov'] >= 0.80 or row['afdb_core_cov'] == 0):
            row['needs_esmf'] = 0
        else:
            row['needs_esmf'] = 1

        # ESMFold fields (placeholder - filled in Phase 3.1A-4)
        row['esmf_plddt_core'] = 0
        row['esmf_core_cov'] = 0

        if afdb_entry:
            time.sleep(args.rate_limit)

    # === Step C: Output panel_candidates.tsv ===
    print("\n=== Step C: Writing panel_candidates.tsv ===")

    output_path = outdir / 'panel_candidates.tsv'
    fields = ['rep_id', 'subtype', 'cluster_id', 'cluster_size', 'seq_len',
              'has_pdb', 'pdb_ids', 'has_afdb', 'afdb_global_plddt',
              'afdb_plddt_core', 'afdb_core_cov',
              'needs_esmf', 'esmf_plddt_core', 'esmf_core_cov']

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r['subtype'], -r['cluster_size'])):
            writer.writerow(row)

    print(f"  Written: {output_path} ({len(rows)} rows)")

    # === Summary ===
    print("\n=== Summary ===")
    for st in ['Ia', 'Ib', 'II']:
        st_rows = [r for r in rows if r['subtype'] == st]
        pdb_count = sum(1 for r in st_rows if r['has_pdb'])
        afdb_count = sum(1 for r in st_rows if r['has_afdb'])
        afdb_good = sum(1 for r in st_rows if r['has_afdb'] and r['afdb_plddt_core'] >= 70)
        needs_esmf = sum(1 for r in st_rows if r['needs_esmf'])
        print(f"  {st}: {len(st_rows)} reps, PDB={pdb_count}, "
              f"AFDB={afdb_count} (core≥70: {afdb_good}), needs_ESMFold={needs_esmf}")

    total_pdb = sum(1 for r in rows if r['has_pdb'])
    total_afdb = sum(1 for r in rows if r['has_afdb'])
    total_esmf = sum(1 for r in rows if r['needs_esmf'])
    print(f"\n  Total: PDB={total_pdb}, AFDB={total_afdb}, needs_ESMFold={total_esmf}")
    print(f"\n  Done. Next: run select_structure_panel.py to validate.")


if __name__ == '__main__':
    main()
