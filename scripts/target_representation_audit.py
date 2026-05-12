#!/usr/bin/env python3
"""Generate Strategy A len255 target representation and rescue summaries."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from statistics import median


TARGETS = [
    ("PfuDAH7PS", "Q8U0A9", "UniRef90_Q8U0A9", 262, "Ib"),
    ("ApeDAH7PS", "Q9YEJ7", "UniRef90_Q9YEJ7", 270, "Ib"),
    ("TmaDAH7PS", "Q9WYH8", "UniRef90_Q9WYH8", 338, "Ib"),
    ("GspDAH7PS", "A0ACD6B8N4", "", 360, "Ib"),
    ("LmoDAH7PS", "Q8Y6T2", "UniRef90_Q8Y6T2", 361, "Ib"),
    ("PniDAH7PS", "V8CS59", "UniRef90_F9DH16", 354, "Ib"),
    ("FtuDAH7PS", "Q5NG89", "UniRef90_Q5NG89", 370, "Ia"),
    ("NmeDAH7PS", "Q9K169", "UniRef90_Q9K169", 351, "Ia"),
    ("SceDAH7PS ARO3", "P14843", "UniRef90_P14843", 370, "Ia"),
    ("SceDAH7PS ARO4", "P32449", "UniRef90_P32449", 370, "Ia"),
    ("EcoDAH7PS AroF", "P00888", "UniRef90_P00888", 356, "Ia"),
    ("EcoDAH7PS AroG", "P0AB91", "UniRef90_P0AB91", 350, "Ia"),
    ("MtuDAH7PS", "O53512", "UniRef90_O53512", 462, "II"),
    ("HpyDAH7PS", "O24947", "UniRef90_O24947", 449, "II"),
    ("CglDAH7PS", "P35170", "UniRef90_P35170", 366, "Ia"),
    ("PaeDAH7PS PA2843", "Q9I000", "UniRef90_Q9I000", 448, "II"),
    ("PaeDAH7PS PA1901", "Q7DC82", "UniRef90_Q7DC82", 405, "II"),
    ("legacy A0A0F2JEB6", "A0A0F2JEB6", "", "", "legacy"),
]

TARGET_COLUMNS = [
    "abbrev",
    "primary_accession",
    "uniref_id",
    "length",
    "subtype_or_expected_group",
    "current_qc_status",
    "len255_qc_status",
    "nr80_status",
    "nr80_representative",
    "nr80_identity_if_absorbed",
    "nr80_local_target_coverage",
    "seeds60_status",
    "final_core_tree_expected_status",
    "label_policy",
    "evidence_source",
    "curation_call",
    "notes",
]


@dataclass(frozen=True)
class ClusterHit:
    cluster: str
    representative: str
    is_representative: bool
    identity: str
    local_target_coverage: str


def read_fasta_headers(path: Path) -> dict[str, str]:
    headers: dict[str, str] = {}
    if not path.exists():
        return headers
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                seq_id = header.split()[0]
                headers[seq_id] = header
    return headers


def read_qc(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {
            row["seq_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("seq_id")
        }


def parse_cdhit_clusters(path: Path) -> dict[str, ClusterHit]:
    if not path.exists():
        return {}
    clusters: list[tuple[str, list[dict[str, object]]]] = []
    current_cluster = ""
    current_members: list[dict[str, object]] = []
    member_re = re.compile(r"^\d+\t(?P<length>\d+)aa, >(?P<seq_id>[^.]+)\.\.\. (?P<tail>.*)$")
    coord_re = re.compile(
        r"at (?P<q_from>\d+):(?P<q_to>\d+):(?P<s_from>\d+):(?P<s_to>\d+)/(?P<identity>[0-9.]+)%"
    )
    simple_identity_re = re.compile(r"at (?P<identity>[0-9.]+)%")

    def flush() -> None:
        if current_cluster:
            clusters.append((current_cluster, list(current_members)))

    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip()
            if line.startswith(">Cluster "):
                flush()
                current_cluster = line.split()[1]
                current_members = []
                continue
            match = member_re.match(line)
            if not match:
                continue
            tail = match.group("tail")
            length = int(match.group("length"))
            coords = coord_re.search(tail)
            coverage = ""
            identity = ""
            if coords:
                q_from = int(coords.group("q_from"))
                q_to = int(coords.group("q_to"))
                aligned = abs(q_to - q_from) + 1
                coverage = f"{aligned / length:.4f}"
                identity = coords.group("identity")
            else:
                simple_identity = simple_identity_re.search(tail)
                if simple_identity:
                    identity = simple_identity.group("identity")
            current_members.append(
                {
                    "seq_id": match.group("seq_id"),
                    "length": length,
                    "tail": tail,
                    "is_rep": tail == "*",
                    "identity": identity,
                    "coverage": coverage,
                }
            )
    flush()

    hits: dict[str, ClusterHit] = {}
    for cluster_id, members in clusters:
        reps = [member["seq_id"] for member in members if member["is_rep"]]
        representative = str(reps[0]) if reps else ""
        for member in members:
            seq_id = str(member["seq_id"])
            is_rep = bool(member["is_rep"])
            hits[seq_id] = ClusterHit(
                cluster=cluster_id,
                representative=representative,
                is_representative=is_rep,
                identity="" if is_rep else str(member["identity"]),
                local_target_coverage="1.0000" if is_rep else str(member["coverage"]),
            )
    return hits


def candidate_ids(primary: str, uniref: str) -> list[str]:
    ids = []
    if uniref:
        ids.append(uniref)
    if primary:
        ids.append(primary)
        if not primary.startswith("UniRef90_"):
            ids.append(f"UniRef90_{primary}")
    return list(dict.fromkeys(ids))


def first_present(ids: list[str], *maps: dict[str, object]) -> str:
    for seq_id in ids:
        if any(seq_id in mapping for mapping in maps):
            return seq_id
    return ""


def status_from_qc(qc: dict[str, dict[str, str]], ids: list[str], missing: str) -> str:
    seq_id = first_present(ids, qc)
    if not seq_id:
        return missing
    row = qc[seq_id]
    return row.get("bin", "present_no_bin")


def classify_header(header: str) -> str:
    text = header.lower()
    dah7ps_terms = [
        "phospho-2-dehydro-3-deoxyheptonate aldolase",
        "2-dehydro-3-deoxyphosphoheptonate aldolase",
        "3-deoxy-7-phosphoheptulonate",
        "3-deoxy-d-arabino-heptulosonate 7-phosphate synthase",
        "3-deoxy-d-arabinoheptulosonate-7-phosphate synthase",
        "3-deoxy-d-arabino-heptulosonate-7-phosphate synthase",
        "d-arabino-heptulosonate-7-phosphate synthase",
        "heptulosonate-7-phosphat e synthase",
        "dah7p",
        "dahp",
        "arof",
    ]
    kdops_terms = [
        "phosphooctonate",
        "manno-octulosonate",
        "kdo8",
        "kdo-8",
    ]
    if any(term in text for term in kdops_terms):
        return "KDOPS-like"
    if "hypothetical" in text or "uncharacterized" in text:
        return "hypothetical/uncharacterized"
    if "chorismate mutase" in text and any(term in text for term in dah7ps_terms):
        return "DAH7PS plus chorismate mutase"
    if any(term in text for term in dah7ps_terms):
        return "DAH7PS-like"
    return "other/ambiguous"


def describe_seed_status(target_id: str, nr80_hit: ClusterHit | None, seed_hits: dict[str, ClusterHit]) -> str:
    if target_id in seed_hits:
        seed_hit = seed_hits[target_id]
        if seed_hit.is_representative:
            return "direct_seeds60_representative"
        return f"seeds60_member_of:{seed_hit.representative};identity={seed_hit.identity}"
    if nr80_hit and nr80_hit.representative in seed_hits:
        seed_hit = seed_hits[nr80_hit.representative]
        if seed_hit.is_representative:
            return f"nr80_representative_is_seeds60_representative:{nr80_hit.representative}"
        return (
            f"nr80_representative_in_seeds60:{nr80_hit.representative}"
            f"->{seed_hit.representative};identity={seed_hit.identity}"
        )
    return "absent_from_seeds60"


def build_target_rows(root: Path) -> list[dict[str, str]]:
    old_qc = {
        "Ib": read_qc(root / "results/02_qc/qc_classification_Ib.tsv"),
        "Ia": read_qc(root / "results/02_qc/qc_classification_Ia.tsv"),
        "II": read_qc(root / "results/02_qc/qc_classification_II.tsv"),
    }
    len255_ib_qc = read_qc(root / "results/02_qc_len255/qc_classification_Ib.tsv")
    nr80_hits = {
        "Ib": parse_cdhit_clusters(root / "results/02_qc_len255/nr80_Ib.fasta.clstr"),
        "Ia": parse_cdhit_clusters(root / "results/02_qc/nr80_Ia.fasta.clstr"),
        "II": parse_cdhit_clusters(root / "results/02_qc/nr80_II.fasta.clstr"),
    }
    seed_hits = {
        "Ib": parse_cdhit_clusters(root / "results/02_qc_len255/seeds60_Ib.fasta.clstr"),
        "Ia": parse_cdhit_clusters(root / "results/02_qc/seeds60_Ia.fasta.clstr"),
        "II": parse_cdhit_clusters(root / "results/02_qc/seeds60_II.fasta.clstr"),
    }
    source_paths = {
        "Ib": (
            "results/02_qc_len255/qc_classification_Ib.tsv;"
            "results/02_qc_len255/nr80_Ib.fasta.clstr;"
            "results/02_qc_len255/seeds60_Ib.fasta.clstr"
        ),
        "Ia": (
            "results/02_qc/qc_classification_Ia.tsv;"
            "results/02_qc/nr80_Ia.fasta.clstr;"
            "results/02_qc/seeds60_Ia.fasta.clstr"
        ),
        "II": (
            "results/02_qc/qc_classification_II.tsv;"
            "results/02_qc/nr80_II.fasta.clstr;"
            "results/02_qc/seeds60_II.fasta.clstr"
        ),
    }

    rows: list[dict[str, str]] = []
    for abbrev, primary, uniref, length, group in TARGETS:
        row = {
            "abbrev": abbrev,
            "primary_accession": primary,
            "uniref_id": uniref,
            "length": str(length),
            "subtype_or_expected_group": group,
            "current_qc_status": "",
            "len255_qc_status": "",
            "nr80_status": "",
            "nr80_representative": "",
            "nr80_identity_if_absorbed": "",
            "nr80_local_target_coverage": "",
            "seeds60_status": "",
            "final_core_tree_expected_status": "",
            "label_policy": "",
            "evidence_source": "",
            "curation_call": "",
            "notes": "",
        }
        if group == "legacy":
            row.update(
                {
                    "current_qc_status": "unresolved_accession",
                    "len255_qc_status": "unresolved_accession",
                    "nr80_status": "unresolved_accession",
                    "seeds60_status": "not_applicable",
                    "final_core_tree_expected_status": "unresolved_accession",
                    "label_policy": "do_not_label_final_tree",
                    "evidence_source": "curation_note",
                    "curation_call": "unresolved_legacy_accession",
                    "notes": "Do not merge silently with V8CS59 / UniRef90_F9DH16.",
                }
            )
            rows.append(row)
            continue

        ids = candidate_ids(primary, uniref)
        qc_for_group = old_qc[group]
        new_qc = len255_ib_qc if group == "Ib" else qc_for_group
        target_id = first_present(ids, new_qc, qc_for_group, nr80_hits[group], seed_hits[group])
        row["current_qc_status"] = status_from_qc(qc_for_group, ids, "absent_from_reference_qc")
        if group == "Ib":
            row["len255_qc_status"] = status_from_qc(new_qc, ids, "absent_from_len255_qc")
        else:
            row["len255_qc_status"] = "not_applicable_non_Ib"
        row["evidence_source"] = source_paths[group]

        if not target_id:
            row.update(
                {
                    "nr80_status": "absent_from_nr80",
                    "seeds60_status": "absent_from_seeds60",
                    "final_core_tree_expected_status": "absent_from_input",
                    "label_policy": "do_not_label_final_tree",
                    "curation_call": "missing_from_staged_reference_artifacts",
                    "notes": "Target identifier not found in staged/formal QC or CD-HIT artifacts.",
                }
            )
            rows.append(row)
            continue

        nr80_hit = nr80_hits[group].get(target_id)
        if nr80_hit is None:
            row.update(
                {
                    "nr80_status": "not_found_in_nr80",
                    "seeds60_status": describe_seed_status(target_id, None, seed_hits[group]),
                    "final_core_tree_expected_status": "needs_curation",
                    "label_policy": "hold_label_pending_curation",
                    "curation_call": "found_before_nr80_but_missing_from_nr80",
                    "notes": f"Resolved target artifact id: {target_id}.",
                }
            )
            rows.append(row)
            continue

        row["nr80_representative"] = nr80_hit.representative
        row["nr80_identity_if_absorbed"] = nr80_hit.identity
        row["nr80_local_target_coverage"] = nr80_hit.local_target_coverage
        row["seeds60_status"] = describe_seed_status(target_id, nr80_hit, seed_hits[group])
        if nr80_hit.is_representative:
            row.update(
                {
                    "nr80_status": "direct_nr80_representative",
                    "final_core_tree_expected_status": "direct_nr80_tip",
                    "label_policy": "label_direct_nr80_tip",
                    "curation_call": "direct_tip_accepted",
                    "notes": f"Resolved target artifact id: {target_id}.",
                }
            )
        else:
            curation_call = "nr80_surrogate_accepted"
            notes = f"Resolved target artifact id: {target_id}; represented by nr80 surrogate."
            if abbrev == "PniDAH7PS":
                curation_call = "nr80_surrogate_tracked_pni_candidate"
                notes = (
                    "Resolved current Pni candidate as UniRef90_F9DH16; "
                    "do not merge with legacy A0A0F2JEB6 without separate curation."
                )
            row.update(
                {
                    "nr80_status": "nr80_cluster_member",
                    "final_core_tree_expected_status": "represented_by_nr80_surrogate",
                    "label_policy": "label_nr80_surrogate_with_target_note",
                    "curation_call": curation_call,
                    "notes": notes,
                }
            )
        rows.append(row)
    return rows


def write_target_table(root: Path, rows: list[dict[str, str]]) -> None:
    out_path = root / "results/02_qc_len255/target_representation.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=TARGET_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def write_rescue_summary(root: Path) -> None:
    old_qc = read_qc(root / "results/02_qc/qc_classification_Ib.tsv")
    new_qc = read_qc(root / "results/02_qc_len255/qc_classification_Ib.tsv")
    headers = read_fasta_headers(root / "results/01_mining/hits_Ib_clean.fasta")
    nr80_headers = read_fasta_headers(root / "results/02_qc_len255/nr80_Ib.fasta")
    seed_headers = read_fasta_headers(root / "results/02_qc_len255/seeds60_Ib.fasta")

    rescued_ids = []
    for seq_id, row in new_qc.items():
        length = int(row["L_seq"])
        if not 255 <= length <= 279 or row.get("bin") != "PASS_CANONICAL":
            continue
        old_bin = old_qc.get(seq_id, {}).get("bin", "missing_old_qc")
        if old_bin == "FRAG":
            rescued_ids.append(seq_id)

    lengths = [int(new_qc[seq_id]["L_seq"]) for seq_id in rescued_ids]
    covs = [float(new_qc[seq_id].get("cov_best", new_qc[seq_id].get("cov_hmm", "0"))) for seq_id in rescued_ids]
    class_counts: dict[str, int] = {}
    class_ids: dict[str, list[str]] = {}
    for seq_id in rescued_ids:
        cls = classify_header(headers.get(seq_id, seq_id))
        class_counts[cls] = class_counts.get(cls, 0) + 1
        class_ids.setdefault(cls, []).append(seq_id)

    def mmx(values: list[int] | list[float], precision: int = 4) -> str:
        if not values:
            return ""
        med = median(values)
        if isinstance(values[0], float):
            return f"min={min(values):.{precision}f};median={med:.{precision}f};max={max(values):.{precision}f}"
        return f"min={min(values)};median={med:g};max={max(values)}"

    metrics = [
        ("255_279_pass_canonical_previously_frag_count", str(len(rescued_ids))),
        ("rescued_entering_nr80_representatives_count", str(sum(seq_id in nr80_headers for seq_id in rescued_ids))),
        ("rescued_entering_seeds60_representatives_count", str(sum(seq_id in seed_headers for seq_id in rescued_ids))),
        ("rescued_length_min_median_max", mmx(lengths, precision=0)),
        ("rescued_cov_best_min_median_max", mmx(covs, precision=4)),
    ]
    for cls in sorted(class_counts):
        metrics.append((f"annotation_class_count:{cls}", str(class_counts[cls])))
    for cls in ["KDOPS-like", "hypothetical/uncharacterized", "other/ambiguous"]:
        ids = class_ids.get(cls, [])
        metrics.append((f"flagged_ids:{cls}", ",".join(ids) if ids else "none"))

    out_path = root / "results/02_qc_len255/len255_rescue_summary.tsv"
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Project root")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    rows = build_target_rows(root)
    write_target_table(root, rows)
    write_rescue_summary(root)


if __name__ == "__main__":
    main()
