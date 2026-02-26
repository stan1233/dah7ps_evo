#!/usr/bin/env python3
"""Augment panel_candidates.tsv with PDB mappings via UniRef90 cluster members.

Phase 3.1A-1 fix: UniRef90 representative may lack PDB, but cluster members
can have PDB crossrefs. This script queries UniProt REST API by cluster ID
and updates has_pdb / pdb_ids accordingly. Optional fallback to RCSB sequence
search is provided for clusters with no PDB crossrefs.

Usage:
  python3 scripts/augment_candidates_with_cluster_pdb.py \
    --candidates_tsv results/03_msa_core/panel_candidates.tsv \
    --out_tsv results/03_msa_core/panel_candidates_clusterpdb.tsv \
    --report_tsv results/03_msa_core/cluster_pdb_report.tsv
"""

import argparse
import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode
import urllib.request
import urllib.error


TRUE_SET = {"1", "true", "yes", "y"}


def parse_bool(value):
    return str(value).strip().lower() in TRUE_SET


def split_ids(value):
    if not value:
        return []
    return [v.strip() for v in str(value).split(',') if v.strip()]


def merge_preserve(a_list, b_list):
    seen = set()
    out = []
    for item in a_list + b_list:
        if not item:
            continue
        item = item.strip()
        if not item:
            continue
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def extract_accession(rep_id):
    acc = rep_id
    if acc.startswith('UniRef90_'):
        acc = acc[len('UniRef90_'):]
    return acc


def is_uniprot_accession(acc):
    return not acc.startswith('UPI')


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def query_uniref_cluster_pdb(
    cluster_id,
    cache_dir,
    max_results=200,
    timeout=20,
    refresh=False,
    retries=2,
    retry_backoff=1.5,
):
    """Query UniProt for members of a UniRef90 cluster with PDB crossrefs."""
    ensure_dir(cache_dir)
    cache_file = cache_dir / f"{cluster_id}_cluster_pdb.json"
    if cache_file.exists() and not refresh:
        with open(cache_file) as f:
            cached = json.load(f)
        # Do not trust transient error cache by default; retry live query.
        if not str(cached.get("status", "")).startswith("error_"):
            return cached

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"uniref_cluster_90:{cluster_id} AND database:pdb",
        "fields": "accession,xref_pdb",
        "format": "json",
        "size": str(max_results),
    }
    url = f"{base_url}?{urlencode(params)}"

    result = {
        "cluster_id": cluster_id,
        "status": "not_found",
        "pdb_ids": [],
        "member_accessions": [],
        "result_count": 0,
        "truncated": False,
        "source": "uniprot_cluster",
    }

    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode())
            results = data.get("results", [])
            pdb_ids = set()
            member_accs = []
            for entry in results:
                acc = entry.get("primaryAccession") or entry.get("accession") or ""
                if acc:
                    member_accs.append(acc)
                xrefs = entry.get("uniProtKBCrossReferences") or entry.get("xrefPdb") or []
                for xref in xrefs:
                    if isinstance(xref, dict) and xref.get("database") == "PDB":
                        pid = xref.get("id")
                        if pid:
                            pdb_ids.add(pid.upper())
            result["pdb_ids"] = sorted(pdb_ids)
            result["member_accessions"] = member_accs
            result["result_count"] = len(results)
            result["truncated"] = len(results) >= int(max_results)
            result["status"] = "found" if pdb_ids else "not_found"
            break
        except urllib.error.HTTPError as e:
            retriable = e.code in {429, 500, 502, 503, 504}
            result["status"] = f"error_{e.code}"
            if retriable and attempt < retries:
                time.sleep(retry_backoff * (2 ** attempt))
                continue
            break
        except Exception as e:
            result["status"] = f"error_{str(e)[:50]}"
            if attempt < retries:
                time.sleep(retry_backoff * (2 ** attempt))
                continue
            break

    with open(cache_file, 'w') as f:
        json.dump(result, f)
    return result


def fetch_uniprot_sequence(accession, cache_dir, timeout=20, refresh=False):
    """Fetch UniProt FASTA and return sequence string."""
    ensure_dir(cache_dir)
    cache_file = cache_dir / f"{accession}.fasta"
    if cache_file.exists() and not refresh:
        seq = []
        with open(cache_file) as f:
            for line in f:
                if line.startswith('>'):
                    continue
                seq.append(line.strip())
        return ''.join(seq)

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        req = urllib.request.Request(url, headers={'Accept': 'text/plain'})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode()
        with open(cache_file, 'w') as f:
            f.write(text)
        seq = []
        for line in text.splitlines():
            if line.startswith('>'):
                continue
            seq.append(line.strip())
        return ''.join(seq)
    except Exception:
        return ""


def query_rcsb_by_sequence(seq, identity_cutoff=0.9, evalue_cutoff=0.1, rows=5, timeout=30):
    """Query RCSB sequence search. Returns polymer_entity identifiers and PDB IDs."""
    if not seq or len(seq) < 25:
        return [], []
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": evalue_cutoff,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": seq,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "paginate": {"start": 0, "rows": rows},
        },
    }
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    data = json.dumps(query).encode('utf-8')
    req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode().strip()
            if not body or getattr(resp, "status", 200) == 204:
                return [], []
            payload = json.loads(body)
    except urllib.error.HTTPError as e:
        if e.code in {204, 404}:
            return [], []
        raise
    except json.JSONDecodeError:
        return [], []
    hits = payload.get("result_set", [])
    entity_ids = []
    pdb_ids = []
    for hit in hits:
        ident = hit.get("identifier")
        if not ident:
            continue
        entity_ids.append(ident)
        pdb_ids.append(ident.split('_')[0].upper())
    return entity_ids, pdb_ids


def load_candidates(path):
    if not Path(path).exists():
        print(f"ERROR: candidates_tsv not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            print("ERROR: Empty candidates TSV", file=sys.stderr)
            sys.exit(1)
        rows = list(reader)
    return rows, list(reader.fieldnames)


def write_tsv(path, rows, fieldnames):
    out_path = Path(path)
    ensure_dir(out_path.parent)
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Augment panel_candidates.tsv with PDB mappings via UniRef90 cluster members."
    )
    parser.add_argument("--candidates_tsv", required=True, help="Input panel_candidates.tsv")
    parser.add_argument("--out_tsv", help="Output TSV (default: add .clusterpdb.tsv)")
    parser.add_argument("--report_tsv", help="Output mapping report TSV")
    parser.add_argument("--cache_dir", help="Cache dir for cluster PDB queries")
    parser.add_argument("--rate_limit", type=float, default=0.2, help="Seconds to sleep between requests")
    parser.add_argument("--max_results", type=int, default=200, help="Max UniProt results per cluster")
    parser.add_argument("--timeout", type=int, default=20, help="HTTP timeout seconds")
    parser.add_argument("--retries", type=int, default=2, help="Retries for transient API errors")
    parser.add_argument("--retry_backoff", type=float, default=1.5, help="Retry backoff base seconds")
    parser.add_argument("--refresh", action="store_true", help="Ignore cache and re-query")
    parser.add_argument("--all", action="store_true", help="Process all rows, not just missing PDB")
    parser.add_argument("--max_rows", type=int, default=0, help="Process first N rows only (0=all)")
    parser.add_argument("--progress_every", type=int, default=25, help="Print progress every N processed rows")
    parser.add_argument("--use_rcsb", action="store_true", help="Fallback to RCSB sequence search")
    parser.add_argument("--rcsb_identity", type=float, default=0.9, help="RCSB identity cutoff (0-1)")
    parser.add_argument("--rcsb_evalue", type=float, default=0.1, help="RCSB evalue cutoff")
    parser.add_argument("--rcsb_rows", type=int, default=5, help="RCSB rows to return")
    args = parser.parse_args()

    rows, fieldnames = load_candidates(args.candidates_tsv)
    extra_fields = ["pdb_member_accs", "pdb_source"]
    for ef in extra_fields:
        if ef not in fieldnames:
            fieldnames.append(ef)

    in_path = Path(args.candidates_tsv)
    out_tsv = args.out_tsv or str(in_path.with_suffix('')) + ".clusterpdb.tsv"
    report_tsv = args.report_tsv or str(in_path.parent / "cluster_pdb_report.tsv")
    cache_dir = Path(args.cache_dir) if args.cache_dir else (in_path.parent / "structure_availability" / "cluster_pdb")

    processed = 0
    added = 0
    report_rows = []
    status_counter = Counter()

    for idx, row in enumerate(rows, start=1):
        if args.max_rows and processed >= args.max_rows:
            break
        rep_id = row.get("rep_id", "")
        cluster_id = row.get("cluster_id", rep_id)
        has_pdb = parse_bool(row.get("has_pdb", "0"))

        if (not args.all) and has_pdb:
            if not row.get("pdb_source"):
                row["pdb_source"] = "sifts" if row.get("pdb_ids") else "none"
            continue

        processed += 1
        mapping = query_uniref_cluster_pdb(
            cluster_id,
            cache_dir,
            max_results=args.max_results,
            timeout=args.timeout,
            refresh=args.refresh,
            retries=args.retries,
            retry_backoff=args.retry_backoff,
        )
        status_counter[mapping.get("status", "unknown")] += 1

        member_accs = mapping.get("member_accessions", [])
        cluster_pdb_ids = mapping.get("pdb_ids", [])

        rcsb_pdb_ids = []
        rcsb_entities = []
        rcsb_used = False
        if args.use_rcsb and not cluster_pdb_ids:
            acc = extract_accession(rep_id)
            if is_uniprot_accession(acc):
                seq_cache = cache_dir / "seq"
                seq = fetch_uniprot_sequence(acc, seq_cache, timeout=args.timeout, refresh=args.refresh)
                if seq:
                    rcsb_entities, rcsb_pdb_ids = query_rcsb_by_sequence(
                        seq,
                        identity_cutoff=args.rcsb_identity,
                        evalue_cutoff=args.rcsb_evalue,
                        rows=args.rcsb_rows,
                        timeout=args.timeout,
                    )
                    rcsb_used = True

        existing_pdb = split_ids(row.get("pdb_ids", ""))
        new_pdb = cluster_pdb_ids if cluster_pdb_ids else rcsb_pdb_ids
        merged_pdb = merge_preserve(existing_pdb, new_pdb)

        if merged_pdb:
            row["has_pdb"] = "1"
            row["pdb_ids"] = ",".join(merged_pdb)
        else:
            row["has_pdb"] = "0"

        if cluster_pdb_ids:
            row["pdb_member_accs"] = ",".join(member_accs)

        # Set pdb_source
        if existing_pdb and cluster_pdb_ids:
            row["pdb_source"] = "sifts+cluster"
        elif existing_pdb:
            row["pdb_source"] = "sifts"
        elif cluster_pdb_ids:
            row["pdb_source"] = "cluster_member"
        elif rcsb_used and rcsb_pdb_ids:
            row["pdb_source"] = "rcsb_seq"
        else:
            row["pdb_source"] = "none"

        if merged_pdb and not existing_pdb:
            added += 1

        report_rows.append({
            "rep_id": rep_id,
            "cluster_id": cluster_id,
            "status": mapping.get("status", ""),
            "pdb_ids": ",".join(cluster_pdb_ids) if cluster_pdb_ids else ",".join(rcsb_pdb_ids),
            "member_accs": ",".join(member_accs),
            "rcsb_entities": ",".join(rcsb_entities),
            "source": row.get("pdb_source", ""),
            "result_count": str(mapping.get("result_count", "")),
            "truncated": str(mapping.get("truncated", "")),
        })

        if args.rate_limit:
            time.sleep(args.rate_limit)
        if args.progress_every and processed % args.progress_every == 0:
            print(
                f"[progress] processed={processed}/{len(rows)} "
                f"added={added} last={rep_id} status={mapping.get('status', '')}",
                flush=True,
            )

    write_tsv(out_tsv, rows, fieldnames)
    if report_rows:
        write_tsv(report_tsv, report_rows, [
            "rep_id", "cluster_id", "status", "pdb_ids", "member_accs",
            "rcsb_entities", "source", "result_count", "truncated",
        ])

    print("=== Cluster PDB augmentation ===")
    print(f"Processed: {processed}")
    print(f"Added PDB mappings: {added}")
    if status_counter:
        print("Status counts:")
        for status, count in sorted(status_counter.items(), key=lambda x: (-x[1], x[0])):
            print(f"  {status}: {count}")
    print(f"Output candidates: {out_tsv}")
    if report_rows:
        print(f"Report: {report_tsv}")


if __name__ == '__main__':
    main()
