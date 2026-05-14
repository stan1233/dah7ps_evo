#!/usr/bin/env python3
"""Round 7A audit for the unrooted Strategy A len255 nr80 core tree.

The script does not run IQ-TREE, root trees, filter alignments, remove
sequences, or run ASR. It reads existing Round 5/6 artifacts and writes only to
results/06_tree_len255/.
"""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/06_tree_len255"

ALIGNMENT = "results/04_msa_len255/core_len255.masked.afa"
ALIGNMENT_QC = "results/04_msa_len255/core_len255.alignment_qc.tsv"
ROUND6_SEQUENCE_CALLS = "results/05_curation_len255/round6_sequence_curation_calls.tsv"
ROUND6_TARGET_REVIEW = "results/05_curation_len255/round6_target_curation_review.tsv"
ROUND6_FLAGGED_REVIEW = "results/05_curation_len255/round6_flagged_rescue_review.tsv"
ROUND6_PNI_REVIEW = "results/05_curation_len255/round6_pni_surrogate_review.tsv"
ROUND6_EXCLUSION_LIST = "results/05_curation_len255/round6_candidate_exclusion_list.tsv"
ROUND6_KEEP_REVIEW = "results/05_curation_len255/round6_candidate_keep_review_list.tsv"
ROUND6_HANDOFF = "results/05_curation_len255/round6_tree_readiness_handoff.md"
SUBTYPE_MAP = "results/02_qc_len255/nr80_all_len255_subtype_map.tsv"
TARGET_REPRESENTATION = "results/02_qc_len255/target_representation.tsv"

TREE_PREFIX = "results/06_tree_len255/core_len255_masked_MFP_unrooted"
TREEFILE = f"{TREE_PREFIX}.treefile"
IQTREE = f"{TREE_PREFIX}.iqtree"
LOG = f"{TREE_PREFIX}.log"
MODEL_GZ = f"{TREE_PREFIX}.model.gz"
CKP_GZ = f"{TREE_PREFIX}.ckp.gz"
COMMAND_FILE = "results/06_tree_len255/round7_iqtree_command.txt"

OUTPUTS = [
    COMMAND_FILE,
    "results/06_tree_len255/round7_tree_summary.tsv",
    "results/06_tree_len255/round7_tipset_check.tsv",
    "results/06_tree_len255/round7_target_tip_presence.tsv",
    "results/06_tree_len255/round7_target_neighborhoods.tsv",
    "results/06_tree_len255/round7_flagged_neighborhoods.tsv",
    "results/06_tree_len255/round7_pni_surrogate_tree_review.tsv",
    "results/06_tree_len255/round7_manual_review_tip_status.tsv",
    "results/06_tree_len255/round7_tree_handoff.md",
    "results/06_tree_len255/round7_artifact_manifest.tsv",
]

IQTREE_OUTPUTS = [TREEFILE, IQTREE, LOG, MODEL_GZ, CKP_GZ]
SOURCE_INPUTS = [
    ALIGNMENT,
    ALIGNMENT_QC,
    "results/04_msa_len255/core_len255_target_carrythrough.tsv",
    "results/04_msa_len255/core_len255_flagged_carrythrough.tsv",
    "results/04_msa_len255/core_len255_subtype_alignment_qc.tsv",
    "results/04_msa_len255/core_len255_artifact_manifest.tsv",
    "results/04_msa_len255/round5_tree_handoff.md",
    "results/05_curation_len255/round6_curation_summary.tsv",
    ROUND6_SEQUENCE_CALLS,
    ROUND6_TARGET_REVIEW,
    ROUND6_FLAGGED_REVIEW,
    ROUND6_PNI_REVIEW,
    ROUND6_EXCLUSION_LIST,
    ROUND6_KEEP_REVIEW,
    ROUND6_HANDOFF,
    TARGET_REPRESENTATION,
    SUBTYPE_MAP,
]

QUERY_TARGETS = [
    ("PfuDAH7PS", "UniRef90_UPI0002AF51CE"),
    ("ApeDAH7PS", "UniRef90_Q9YEJ7"),
    ("PniDAH7PS", "UniRef90_A0A379DXQ3"),
]
FLAGGED_QUERY = "UniRef90_A0A5E4LPI9"


def read_tsv(relative_path: str) -> list[dict[str, str]]:
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(relative_path: str, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path = ROOT / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean_value(row.get(field, "")) for field in fieldnames})


def clean_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def read_text(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


def raw_md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_bool(args: list[str], relative_path: str) -> bool:
    result = subprocess.run(
        ["git", *args, relative_path],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def is_tracked(relative_path: str) -> bool:
    return git_bool(["ls-files", "--error-unmatch"], relative_path)


def is_ignored(relative_path: str) -> bool:
    return git_bool(["check-ignore", "-q"], relative_path)


def alignment_tips(relative_path: str) -> list[str]:
    tips: list[str] = []
    with (ROOT / relative_path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                tips.append(line[1:].strip().split()[0])
    return tips


def alignment_columns(relative_path: str) -> int:
    current: list[str] = []
    lengths: list[int] = []
    with (ROOT / relative_path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if current:
                    lengths.append(len("".join(current)))
                    current = []
            else:
                current.append(line.strip())
        if current:
            lengths.append(len("".join(current)))
    unique_lengths = set(lengths)
    if len(unique_lengths) != 1:
        raise ValueError(f"non-uniform alignment lengths: {sorted(unique_lengths)[:5]}")
    return lengths[0]


class ParsedTree:
    def __init__(self) -> None:
        self.adj: dict[int, list[tuple[int, float]]] = defaultdict(list)
        self.tip_to_node: dict[str, int] = {}
        self.node_to_tip: dict[int, str] = {}
        self.next_id = 0

    def new_node(self, tip_label: str | None = None) -> int:
        node_id = self.next_id
        self.next_id += 1
        if tip_label:
            self.tip_to_node[tip_label] = node_id
            self.node_to_tip[node_id] = tip_label
        return node_id

    def add_edge(self, parent: int, child: int, length: float) -> None:
        self.adj[parent].append((child, length))
        self.adj[child].append((parent, length))


def parse_newick(relative_path: str) -> ParsedTree:
    text = (ROOT / relative_path).read_text(encoding="utf-8").strip()
    parser = NewickParser(text)
    return parser.parse()


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.i = 0
        self.tree = ParsedTree()

    def parse(self) -> ParsedTree:
        self._parse_subtree()
        self._skip_ws()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        self._skip_ws()
        if self.i != len(self.text):
            raise ValueError(f"unexpected Newick suffix at offset {self.i}")
        return self.tree

    def _parse_subtree(self) -> tuple[int, float]:
        self._skip_ws()
        if self.text[self.i] == "(":
            self.i += 1
            node_id = self.tree.new_node()
            while True:
                child_id, child_length = self._parse_subtree()
                self.tree.add_edge(node_id, child_id, child_length)
                self._skip_ws()
                if self.text[self.i] == ",":
                    self.i += 1
                    continue
                if self.text[self.i] == ")":
                    self.i += 1
                    break
                raise ValueError(f"unexpected Newick char {self.text[self.i]!r} at {self.i}")
            self._parse_label()
            length = self._parse_length()
            return node_id, length
        label = self._parse_label()
        if not label:
            raise ValueError(f"empty leaf label at offset {self.i}")
        node_id = self.tree.new_node(label)
        length = self._parse_length()
        return node_id, length

    def _skip_ws(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _parse_label(self) -> str:
        self._skip_ws()
        if self.i < len(self.text) and self.text[self.i] == "'":
            self.i += 1
            start = self.i
            while self.i < len(self.text) and self.text[self.i] != "'":
                self.i += 1
            label = self.text[start : self.i]
            if self.i < len(self.text) and self.text[self.i] == "'":
                self.i += 1
            return label
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,();":
            self.i += 1
        return self.text[start : self.i].strip()

    def _parse_length(self) -> float:
        self._skip_ws()
        if self.i >= len(self.text) or self.text[self.i] != ":":
            return 0.0
        self.i += 1
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",();":
            self.i += 1
        raw = self.text[start : self.i].strip()
        return float(raw) if raw else 0.0


def nearest_neighbors(tree: ParsedTree, query_tip: str, limit: int = 20) -> list[tuple[str, float]]:
    if query_tip not in tree.tip_to_node:
        return []
    start = tree.tip_to_node[query_tip]
    stack = [(start, -1, 0.0)]
    neighbors: list[tuple[str, float]] = []
    while stack:
        node_id, parent, distance = stack.pop()
        tip = tree.node_to_tip.get(node_id)
        if tip and tip != query_tip:
            neighbors.append((tip, distance))
        for next_id, branch_length in tree.adj[node_id]:
            if next_id != parent:
                stack.append((next_id, node_id, distance + branch_length))
    return sorted(neighbors, key=lambda item: (item[1], item[0]))[:limit]


def model_selected() -> str:
    path = ROOT / IQTREE
    if not path.is_file():
        return ""
    text = path.read_text(encoding="utf-8", errors="replace")
    for pattern in [
        r"Best-fit model according to BIC:\s*(\S+)",
        r"Model of substitution:\s*(\S+)",
    ]:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
    return ""


def iqtree_completed() -> bool:
    required = [ROOT / TREEFILE, ROOT / IQTREE, ROOT / LOG]
    if not all(path.is_file() and path.stat().st_size > 0 for path in required):
        return False
    text = read_text(LOG)
    markers = ["Total CPU time used", "Date and Time", "Analysis results written to"]
    return any(marker in text for marker in markers)


def load_maps() -> dict[str, object]:
    sequence_rows = read_tsv(ROUND6_SEQUENCE_CALLS)
    target_rows = read_tsv(ROUND6_TARGET_REVIEW)
    flagged_rows = read_tsv(ROUND6_FLAGGED_REVIEW)
    pni_rows = read_tsv(ROUND6_PNI_REVIEW)
    subtype_rows = read_tsv(SUBTYPE_MAP)
    alignment_qc_rows = read_tsv(ALIGNMENT_QC)
    return {
        "sequence_by_id": {row["seq_id"]: row for row in sequence_rows},
        "target_rows": target_rows,
        "flagged_rows": flagged_rows,
        "pni_rows": pni_rows,
        "subtype_by_id": {row["seq_id"]: row["subtype_source"] for row in subtype_rows},
        "alignment_qc_by_id": {row["seq_id"]: row for row in alignment_qc_rows},
    }


def context_for_tip(seq_id: str, maps: dict[str, object]) -> dict[str, str]:
    sequence_by_id = maps["sequence_by_id"]
    subtype_by_id = maps["subtype_by_id"]
    assert isinstance(sequence_by_id, dict)
    assert isinstance(subtype_by_id, dict)
    row = sequence_by_id.get(seq_id, {})
    return {
        "subtype": subtype_by_id.get(seq_id, row.get("subtype_source", "unknown")),
        "round6_call": row.get("round6_curation_call", ""),
        "is_target": row.get("is_target_or_representative", "false"),
        "is_flagged": row.get("is_flagged_rescue", "false"),
    }


def build_tipset_check(aln_tips: list[str], tree_tips: list[str]) -> list[dict[str, object]]:
    aln_counts = Counter(aln_tips)
    tree_counts = Counter(tree_tips)
    aln_set = set(aln_tips)
    tree_set = set(tree_tips)
    missing = sorted(aln_set - tree_set)
    extra = sorted(tree_set - aln_set)
    duplicate_aln = sorted([tip for tip, count in aln_counts.items() if count > 1])
    duplicate_tree = sorted([tip for tip, count in tree_counts.items() if count > 1])
    kdops_tips = sorted([tip for tip in tree_tips if "KDOPS" in tip])
    return [
        {"metric": "alignment_tip_count", "value": len(aln_tips), "notes": "parsed from masked AFA"},
        {"metric": "tree_tip_count", "value": len(tree_tips), "notes": "parsed from unrooted treefile"},
        {
            "metric": "tip_count_matches",
            "value": len(aln_tips) == len(tree_tips),
            "notes": "count equality only",
        },
        {
            "metric": "duplicate_alignment_ids_count",
            "value": len(duplicate_aln),
            "notes": ",".join(duplicate_aln),
        },
        {
            "metric": "duplicate_tree_ids_count",
            "value": len(duplicate_tree),
            "notes": ",".join(duplicate_tree),
        },
        {"metric": "missing_from_tree_count", "value": len(missing), "notes": "alignment tips absent from tree"},
        {"metric": "extra_in_tree_count", "value": len(extra), "notes": "tree tips absent from alignment"},
        {"metric": "missing_from_tree_ids", "value": ",".join(missing), "notes": "empty means none"},
        {"metric": "extra_in_tree_ids", "value": ",".join(extra), "notes": "empty means none"},
        {
            "metric": "KDOPS_O66496_present",
            "value": "KDOPS_O66496" in tree_set,
            "notes": "expected false for ingroup nr80 len255 tree",
        },
        {
            "metric": "any_KDOPS_tip_present",
            "value": bool(kdops_tips),
            "notes": ",".join(kdops_tips),
        },
    ]


def label_policy_for(row: dict[str, str]) -> str:
    label = row["target_label"]
    if label == "ApeDAH7PS":
        return "ApeDAH7PS direct nr80 tip"
    if label == "PfuDAH7PS":
        return "PfuDAH7PS-like / represents Q8U0A9"
    return row["label_policy"]


def target_tree_status(row: dict[str, str], tree_tip_set: set[str]) -> str:
    label = row["target_label"]
    expected_tip = row["alignment_seq_id"] or row["representative_id"]
    present = expected_tip in tree_tip_set if expected_tip else False
    if label == "legacy A0A0F2JEB6":
        return "unresolved_accession_not_used"
    if label == "PniDAH7PS":
        return "present_but_needs_curation" if present else "missing_needs_curation"
    if row["expected_status_from_target_representation"] == "direct_nr80_tip":
        if present and row["curation_call"] == "manual_review_before_tree":
            return "direct_nr80_tip_present_manual_review_flag_active"
        return "direct_nr80_tip_present" if present else "direct_nr80_tip_missing"
    if row["expected_status_from_target_representation"] == "represented_by_nr80_surrogate":
        return "nr80_surrogate_present" if present else "nr80_surrogate_missing"
    return row["expected_status_from_target_representation"]


def build_target_presence(target_rows: list[dict[str, str]], tree_tip_set: set[str]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for row in target_rows:
        expected_tip = row["alignment_seq_id"] or row["representative_id"]
        expected_status = row["expected_status_from_target_representation"]
        if row["target_label"] == "PniDAH7PS":
            expected_status = "represented_by_nr80_surrogate + needs_curation"
        output.append(
            {
                "target_label": row["target_label"],
                "primary_accession": row["primary_accession"],
                "expected_status_from_target_representation": expected_status,
                "representative_id": row["representative_id"],
                "expected_tree_tip_id": expected_tip,
                "present_in_alignment": row["present_in_alignment"],
                "present_in_tree": expected_tip in tree_tip_set if expected_tip else False,
                "round6_curation_call": row["curation_call"],
                "round7_tree_status": target_tree_status(row, tree_tip_set),
                "label_policy": label_policy_for(row),
                "notes": row["notes"],
            }
        )
    return output


def neighborhood_rows(
    tree: ParsedTree,
    maps: dict[str, object],
    query_label: str,
    query_tip_id: str,
    flagged: bool = False,
    flag_class: str = "",
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rank, (neighbor, distance) in enumerate(nearest_neighbors(tree, query_tip_id), start=1):
        ctx = context_for_tip(neighbor, maps)
        base = {
            "query_label": query_label,
            "query_tip_id": query_tip_id,
            "rank": rank,
            "neighbor_tip_id": neighbor,
            "neighbor_subtype_source": ctx["subtype"],
            "patristic_distance": f"{distance:.8f}",
            "neighbor_round6_curation_call": ctx["round6_call"],
            "neighbor_is_target_or_representative": ctx["is_target"],
            "neighbor_is_flagged_rescue": ctx["is_flagged"],
            "notes": "nearest-neighbor context in the unrooted Round 7A tree",
        }
        if flagged:
            base.update(
                {
                    "flag_class": flag_class,
                    "round7_flagged_status": "present_manual_review_before_final_tree",
                    "recommended_next_action": "manual_review_before_final_tree",
                }
            )
        rows.append(base)
    return rows


def subtype_summary(rows: list[dict[str, object]]) -> str:
    counts = Counter(str(row["neighbor_subtype_source"]) for row in rows)
    return ";".join(f"{key}:{counts[key]}" for key in sorted(counts))


def build_pni_review(
    pni_rows: list[dict[str, str]],
    target_neighborhoods: list[dict[str, object]],
    tree_tip_set: set[str],
) -> list[dict[str, object]]:
    pni_by_label = {row["candidate_label"]: row for row in pni_rows}
    rep_tip = "UniRef90_A0A379DXQ3"
    pni_neighbor_rows = [row for row in target_neighborhoods if row["query_tip_id"] == rep_tip]
    neighbor_summary = subtype_summary(pni_neighbor_rows)
    base = pni_by_label.get("V8CS59 / UniRef90_F9DH16", {})
    rep = pni_by_label.get("UniRef90_A0A379DXQ3", {})
    legacy = pni_by_label.get("legacy A0A0F2JEB6", {})
    return [
        {
            "candidate_label": "V8CS59 / UniRef90_F9DH16",
            "primary_accession": "V8CS59",
            "representative_id": rep_tip,
            "present_in_alignment": base.get("present_in_alignment", ""),
            "present_in_tree": rep_tip in tree_tip_set,
            "gap_fraction": base.get("gap_fraction", ""),
            "round6_curation_call": base.get("curation_call", "needs_curation"),
            "nearest_neighbor_subtypes_summary": neighbor_summary,
            "annotation_status": base.get("annotation_status", "needs_curation"),
            "taxonomy_status": base.get("taxonomy_status", "需要进一步验证"),
            "coverage_status": base.get("coverage_status", "需要进一步验证"),
            "round7_tree_interpretation": "present_but_tree_presence_does_not_resolve_pni_acceptability",
            "recommended_next_action": "keep_needs_curation_pending_annotation_taxonomy_review",
            "notes": "current Pni candidate remains represented by nr80 surrogate; not accepted in Round 7A",
        },
        {
            "candidate_label": "UniRef90_A0A379DXQ3",
            "primary_accession": "A0A379DXQ3",
            "representative_id": rep_tip,
            "present_in_alignment": rep.get("present_in_alignment", ""),
            "present_in_tree": rep_tip in tree_tip_set,
            "gap_fraction": rep.get("gap_fraction", ""),
            "round6_curation_call": rep.get("curation_call", "needs_curation"),
            "nearest_neighbor_subtypes_summary": neighbor_summary,
            "annotation_status": rep.get("annotation_status", "representative_for_Pni_candidate_needs_curation"),
            "taxonomy_status": rep.get("taxonomy_status", "需要进一步验证"),
            "coverage_status": rep.get("coverage_status", "需要进一步验证"),
            "round7_tree_interpretation": "representative_present_but_needs_curation",
            "recommended_next_action": "do_not_accept_surrogate_without_formal_curation",
            "notes": "nearest-neighbor context is unrooted and representation-focused",
        },
        {
            "candidate_label": "legacy A0A0F2JEB6",
            "primary_accession": "A0A0F2JEB6",
            "representative_id": legacy.get("representative_id", ""),
            "present_in_alignment": legacy.get("present_in_alignment", "false"),
            "present_in_tree": False,
            "gap_fraction": legacy.get("gap_fraction", ""),
            "round6_curation_call": legacy.get("curation_call", "needs_external_validation"),
            "nearest_neighbor_subtypes_summary": "",
            "annotation_status": "unresolved_accession",
            "taxonomy_status": legacy.get("taxonomy_status", "需要进一步验证"),
            "coverage_status": legacy.get("coverage_status", "需要进一步验证"),
            "round7_tree_interpretation": "unresolved_accession_not_used",
            "recommended_next_action": "requires_external_validation_if_reintroduced",
            "notes": "not merged with V8CS59 / UniRef90_F9DH16",
        },
    ]


def build_manual_review_rows(sequence_by_id: dict[str, dict[str, str]], tree_tip_set: set[str]) -> list[dict[str, object]]:
    rows = []
    for seq_id in sorted(sequence_by_id):
        row = sequence_by_id[seq_id]
        if row["round6_curation_call"] != "manual_review_before_tree":
            continue
        present = seq_id in tree_tip_set
        rows.append(
            {
                "seq_id": seq_id,
                "subtype_source": row["subtype_source"],
                "gap_fraction": row["gap_fraction"],
                "gap_or_ambiguity_fraction": row["gap_or_ambiguity_fraction"],
                "present_in_tree": present,
                "is_target_or_representative": row["is_target_or_representative"],
                "target_label_if_any": row["target_label_if_any"],
                "is_flagged_rescue": row["is_flagged_rescue"],
                "flag_class_if_any": row["flag_class_if_any"],
                "round6_reason": row["reason"],
                "round7_status": "present_as_round6_manual_review_flag" if present else "missing_from_tree",
                "notes": "Round 6 manual-review call remains an active interpretation flag.",
            }
        )
    return rows


def value_for_metric(rows: list[dict[str, object]], metric: str) -> str:
    for row in rows:
        if row["metric"] == metric:
            return str(row["value"])
    return ""


def metric_is_true(rows: list[dict[str, object]], metric: str) -> bool:
    return value_for_metric(rows, metric).lower() == "true"


def metric_int(rows: list[dict[str, object]], metric: str) -> int:
    raw = value_for_metric(rows, metric)
    return int(raw) if raw else 0


def tipset_is_exact(rows: list[dict[str, object]]) -> bool:
    return (
        metric_int(rows, "alignment_tip_count") == 9677
        and metric_int(rows, "tree_tip_count") == 9677
        and metric_int(rows, "duplicate_alignment_ids_count") == 0
        and metric_int(rows, "duplicate_tree_ids_count") == 0
        and metric_int(rows, "missing_from_tree_count") == 0
        and metric_int(rows, "extra_in_tree_count") == 0
    )


def target_status(rows: list[dict[str, object]], label: str) -> str:
    for row in rows:
        if row["target_label"] == label:
            return str(row["round7_tree_status"])
    return "missing_from_round7_target_presence"


def flagged_status(rows: list[dict[str, object]], seq_id: str) -> str:
    for row in rows:
        if row["query_tip_id"] == seq_id:
            return "present_manual_review_before_final_tree"
    return "missing_from_flagged_neighborhoods"


def build_summary(
    aln_count: int,
    aln_columns: int,
    tipset_rows: list[dict[str, object]],
    target_rows: list[dict[str, object]],
    flagged_rows: list[dict[str, object]],
    manual_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    command = read_text(COMMAND_FILE).strip()
    return [
        {"metric": "alignment_used", "value": ALIGNMENT, "notes": "Round 5 masked alignment"},
        {"metric": "alignment_is_masked", "value": "true", "notes": "RF-masked AFA"},
        {"metric": "alignment_sequence_count", "value": aln_count, "notes": "alignment tip count"},
        {"metric": "alignment_columns", "value": aln_columns, "notes": "alignment column count"},
        {"metric": "tree_prefix", "value": TREE_PREFIX, "notes": "unrooted IQ-TREE prefix"},
        {"metric": "iqtree_binary", "value": "conda run -n dah7ps_v4 iqtree", "notes": "IQ-TREE 3.0.1"},
        {"metric": "iqtree_command", "value": command, "notes": "no -o and no support flags"},
        {"metric": "iqtree_completed", "value": iqtree_completed(), "notes": "treefile/iqtree/log completion check"},
        {"metric": "model_selected", "value": model_selected(), "notes": "best-fit model from .iqtree"},
        {"metric": "tree_tip_count", "value": value_for_metric(tipset_rows, "tree_tip_count"), "notes": "parsed from treefile"},
        {
            "metric": "tipset_matches_alignment",
            "value": value_for_metric(tipset_rows, "missing_from_tree_count") == "0"
            and value_for_metric(tipset_rows, "extra_in_tree_count") == "0",
            "notes": "exact set equality",
        },
        {
            "metric": "KDOPS_O66496_present",
            "value": value_for_metric(tipset_rows, "KDOPS_O66496_present"),
            "notes": "expected false",
        },
        {
            "metric": "any_KDOPS_tip_present",
            "value": value_for_metric(tipset_rows, "any_KDOPS_tip_present"),
            "notes": "expected false",
        },
        {
            "metric": "Q8U0A9_representative_tree_status",
            "value": target_status(target_rows, "PfuDAH7PS"),
            "notes": "Q8U0A9 represented by UniRef90_UPI0002AF51CE",
        },
        {
            "metric": "Q9YEJ7_tree_status",
            "value": target_status(target_rows, "ApeDAH7PS"),
            "notes": "direct nr80 tip with Round 6 manual-review flag",
        },
        {
            "metric": "Pni_candidate_representative_tree_status",
            "value": target_status(target_rows, "PniDAH7PS"),
            "notes": "Pni surrogate remains needs_curation",
        },
        {
            "metric": "legacy_A0A0F2JEB6_tree_status",
            "value": target_status(target_rows, "legacy A0A0F2JEB6"),
            "notes": "unresolved accession not used",
        },
        {
            "metric": "KDOPS_like_UniRef90_A0A5E4LPI9_tree_status",
            "value": flagged_status(flagged_rows, FLAGGED_QUERY),
            "notes": "flagged hit remains manual_review_before_final_tree",
        },
        {
            "metric": "manual_review_before_tree_tip_count",
            "value": sum(row["present_in_tree"] == "true" or row["present_in_tree"] is True for row in manual_rows),
            "notes": "Round 6 manual-review tips present in tree",
        },
        {
            "metric": "candidate_exclude_before_tree_applied_count",
            "value": 0,
            "notes": "Round 6 candidate exclusions were advisory only",
        },
        {"metric": "sequence_removals_applied", "value": "false", "notes": "alignment used as-is"},
        {"metric": "filtered_alignment_created", "value": "false", "notes": "no filtered alignment written"},
        {
            "metric": "round7_interpretation_scope",
            "value": "unrooted_nr80_len255_representation_tree_only",
            "notes": "no rooting, QC3, ASR, or root-sensitive claim",
        },
        {
            "metric": "next_recommended_action",
            "value": "user_review_unrooted_placement_audits_before_any_rooting_or_QC3_decision",
            "notes": "Round 6 manual-review calls remain active",
        },
    ]


def render_handoff(
    summary_rows: list[dict[str, object]],
    target_rows: list[dict[str, object]],
    target_neighborhoods: list[dict[str, object]],
    manual_rows: list[dict[str, object]],
    flagged_rows: list[dict[str, object]],
) -> str:
    summary = {str(row["metric"]): str(row["value"]) for row in summary_rows}
    target_by_label = {str(row["target_label"]): row for row in target_rows}
    query_summaries = []
    for label, tip in QUERY_TARGETS:
        rows = [row for row in target_neighborhoods if row["query_tip_id"] == tip]
        query_summaries.append(f"- {label} / `{tip}`: {subtype_summary(rows)}")
    flagged_status_text = flagged_status(flagged_rows, FLAGGED_QUERY)
    return f"""# Strategy A Round 7A Tree Handoff

This is an unrooted Strategy A len255 nr80 representation tree.
No rooting, QC3, ASR, or root-sensitive evolutionary claim is released by this round.
Round 6 manual-review calls remain active interpretation flags.
No sequence exclusions were applied.

## 1. Tree Run

- Alignment: `{ALIGNMENT}`
- Masking: masked
- Sequences: {summary['alignment_sequence_count']}
- Columns: {summary['alignment_columns']}
- IQ-TREE command: `{summary['iqtree_command']}`
- Completed cleanly: {summary['iqtree_completed']}
- Model selected: {summary['model_selected']}
- Tree prefix: `{TREE_PREFIX}`

## 2. Tip Set

- Tree tip count: {summary['tree_tip_count']}
- Tip set matches alignment: {summary['tipset_matches_alignment']}
- Any KDOPS tip present: {summary['any_KDOPS_tip_present']}
- KDOPS_O66496 present: {summary['KDOPS_O66496_present']}

In the original MFP/LGC20 O66496-containing trees, KDOPS_O66496 is nearest to
Ib-labelled ingroup tips.

## 3. Target Status

- Q8U0A9 / PfuDAH7PS representation: {target_by_label['PfuDAH7PS']['round7_tree_status']} via `UniRef90_UPI0002AF51CE`
- Q9YEJ7 / ApeDAH7PS: {target_by_label['ApeDAH7PS']['round7_tree_status']} via `UniRef90_Q9YEJ7`
- V8CS59 / UniRef90_F9DH16 / PniDAH7PS candidate: {target_by_label['PniDAH7PS']['round7_tree_status']} via `UniRef90_A0A379DXQ3`
- legacy A0A0F2JEB6: {target_by_label['legacy A0A0F2JEB6']['round7_tree_status']}
- KDOPS-like `UniRef90_A0A5E4LPI9`: {flagged_status_text}

The Pni surrogate remains `needs_curation`; tree presence does not resolve
annotation, taxonomy, or acceptability.

## 4. Nearest-Neighbor Context

Nearest-neighbor audits report patristic neighbors in the unrooted Round 7A tree:

{chr(10).join(query_summaries)}

Full tables:

- `results/06_tree_len255/round7_target_neighborhoods.tsv`
- `results/06_tree_len255/round7_flagged_neighborhoods.tsv`
- `results/06_tree_len255/round7_pni_surrogate_tree_review.tsv`

## 5. Manual-Review Tips

Round 6 manual-review tips present in the tree: {summary['manual_review_before_tree_tip_count']}.
These calls remain interpretation flags, not exclusion rules.

## 6. Interpretation Scope

This tree can support unrooted nr80 len255 representation and local
nearest-neighbor placement context for the audited targets and flagged hits.
It does not support rooting, root stability, root-sensitive ASR, or final
directional evolutionary interpretation.

## 7. Remains HOLD

- noO66496 formal S1
- formal S2 noO66496
- QC3 root stability
- root-sensitive ASR
- Pni surrogate acceptance
- KDOPS-like flagged hit acceptance

## 8. Must Not Be Claimed

- solved rooting
- released QC3/root-stability status
- released ASR status
- direct-tip status for every target
- accepted Pni surrogate
- Type Ia placement for O66496
"""


def build_manifest() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    source_inputs = ";".join(SOURCE_INPUTS)
    paths: list[tuple[str, str]] = []
    paths.append((ALIGNMENT, "round5_masked_alignment_input"))
    paths.extend((path, "iqtree_output") for path in IQTREE_OUTPUTS)
    paths.extend((path, "round7_output") for path in OUTPUTS)
    paths.append(("scripts/round7_len255_tree_audit.py", "round7_helper_script"))
    for relative_path in SOURCE_INPUTS:
        if relative_path != ALIGNMENT:
            paths.append((relative_path, "round5_round6_source_input"))

    seen: set[str] = set()
    for relative_path, role in paths:
        if relative_path in seen:
            continue
        seen.add(relative_path)
        path = ROOT / relative_path
        if path.exists() and relative_path != "results/06_tree_len255/round7_artifact_manifest.tsv":
            size = path.stat().st_size
            md5 = raw_md5_file(path)
        elif relative_path == "results/06_tree_len255/round7_artifact_manifest.tsv":
            size = path.stat().st_size if path.exists() else ""
            md5 = "SELF_REFERENTIAL"
        else:
            size = ""
            md5 = "MISSING"
        if role == "iqtree_output":
            required = relative_path in {TREEFILE, IQTREE, LOG}
            recreate = read_text(COMMAND_FILE).strip()
            notes = "IQ-TREE output family; large artifacts may be ignored"
        elif role in {"round7_output", "round7_helper_script"}:
            required = True
            recreate = "python scripts/round7_len255_tree_audit.py"
            notes = "Round 7A audit output or helper"
        else:
            required = True
            recreate = "source artifact from Round 5 or Round 6"
            notes = "Input used by Round 7A"
        if relative_path == "results/06_tree_len255/round7_artifact_manifest.tsv":
            notes = "Self-referential manifest row; md5 cannot include its own final md5"
        rows.append(
            {
                "relative_path": relative_path,
                "artifact_role": role,
                "file_size_bytes": size,
                "md5": md5,
                "tracked_by_git": is_tracked(relative_path),
                "ignored_by_git": is_ignored(relative_path),
                "required_for_next_stage": required,
                "source_inputs": source_inputs if role in {"iqtree_output", "round7_output", "round7_helper_script"} else "",
                "recreate_command_or_source": recreate,
                "notes": notes,
            }
        )
    return rows


def main() -> int:
    for relative_path in [ALIGNMENT, TREEFILE, IQTREE, LOG, COMMAND_FILE, *SOURCE_INPUTS]:
        if not (ROOT / relative_path).is_file():
            raise FileNotFoundError(relative_path)
    if not iqtree_completed():
        raise RuntimeError("IQ-TREE output is incomplete; refusing to write completed Round 7A audit")

    maps = load_maps()
    sequence_by_id = maps["sequence_by_id"]
    target_rows = maps["target_rows"]
    flagged_review_rows = maps["flagged_rows"]
    pni_review_rows = maps["pni_rows"]
    assert isinstance(sequence_by_id, dict)
    assert isinstance(target_rows, list)
    assert isinstance(flagged_review_rows, list)
    assert isinstance(pni_review_rows, list)

    aln_tips = alignment_tips(ALIGNMENT)
    aln_columns = alignment_columns(ALIGNMENT)
    tree = parse_newick(TREEFILE)
    tree_tips = sorted(tree.tip_to_node)
    tree_tip_set = set(tree_tips)
    tipset_rows = build_tipset_check(aln_tips, tree_tips)

    if metric_is_true(tipset_rows, "any_KDOPS_tip_present"):
        write_tsv(
            "results/06_tree_len255/round7_tipset_check.tsv",
            ["metric", "value", "notes"],
            tipset_rows,
        )
        raise RuntimeError("KDOPS tip present in ingroup nr80 tree; stop before placement audit")
    if not tipset_is_exact(tipset_rows):
        write_tsv(
            "results/06_tree_len255/round7_tipset_check.tsv",
            ["metric", "value", "notes"],
            tipset_rows,
        )
        raise RuntimeError("Tree tip set does not exactly match the 9677-tip alignment")

    target_presence_rows = build_target_presence(target_rows, tree_tip_set)
    target_neighborhood_rows: list[dict[str, object]] = []
    for query_label, query_tip_id in QUERY_TARGETS:
        target_neighborhood_rows.extend(neighborhood_rows(tree, maps, query_label, query_tip_id))

    flag_by_seq = {row["seq_id"]: row for row in flagged_review_rows}
    flagged_neighborhood_rows = neighborhood_rows(
        tree,
        maps,
        "KDOPS-like flagged rescue hit",
        FLAGGED_QUERY,
        flagged=True,
        flag_class=flag_by_seq.get(FLAGGED_QUERY, {}).get("flag_class", "KDOPS-like"),
    )
    pni_tree_review_rows = build_pni_review(pni_review_rows, target_neighborhood_rows, tree_tip_set)
    manual_review_rows = build_manual_review_rows(sequence_by_id, tree_tip_set)
    summary_rows = build_summary(
        len(aln_tips),
        aln_columns,
        tipset_rows,
        target_presence_rows,
        flagged_neighborhood_rows,
        manual_review_rows,
    )

    write_tsv("results/06_tree_len255/round7_tipset_check.tsv", ["metric", "value", "notes"], tipset_rows)
    write_tsv(
        "results/06_tree_len255/round7_target_tip_presence.tsv",
        [
            "target_label",
            "primary_accession",
            "expected_status_from_target_representation",
            "representative_id",
            "expected_tree_tip_id",
            "present_in_alignment",
            "present_in_tree",
            "round6_curation_call",
            "round7_tree_status",
            "label_policy",
            "notes",
        ],
        target_presence_rows,
    )
    write_tsv(
        "results/06_tree_len255/round7_target_neighborhoods.tsv",
        [
            "query_label",
            "query_tip_id",
            "rank",
            "neighbor_tip_id",
            "neighbor_subtype_source",
            "patristic_distance",
            "neighbor_round6_curation_call",
            "neighbor_is_target_or_representative",
            "neighbor_is_flagged_rescue",
            "notes",
        ],
        target_neighborhood_rows,
    )
    write_tsv(
        "results/06_tree_len255/round7_flagged_neighborhoods.tsv",
        [
            "query_label",
            "query_tip_id",
            "rank",
            "neighbor_tip_id",
            "neighbor_subtype_source",
            "patristic_distance",
            "neighbor_round6_curation_call",
            "neighbor_is_target_or_representative",
            "neighbor_is_flagged_rescue",
            "notes",
            "flag_class",
            "round7_flagged_status",
            "recommended_next_action",
        ],
        flagged_neighborhood_rows,
    )
    write_tsv(
        "results/06_tree_len255/round7_pni_surrogate_tree_review.tsv",
        [
            "candidate_label",
            "primary_accession",
            "representative_id",
            "present_in_alignment",
            "present_in_tree",
            "gap_fraction",
            "round6_curation_call",
            "nearest_neighbor_subtypes_summary",
            "annotation_status",
            "taxonomy_status",
            "coverage_status",
            "round7_tree_interpretation",
            "recommended_next_action",
            "notes",
        ],
        pni_tree_review_rows,
    )
    write_tsv(
        "results/06_tree_len255/round7_manual_review_tip_status.tsv",
        [
            "seq_id",
            "subtype_source",
            "gap_fraction",
            "gap_or_ambiguity_fraction",
            "present_in_tree",
            "is_target_or_representative",
            "target_label_if_any",
            "is_flagged_rescue",
            "flag_class_if_any",
            "round6_reason",
            "round7_status",
            "notes",
        ],
        manual_review_rows,
    )
    write_tsv("results/06_tree_len255/round7_tree_summary.tsv", ["metric", "value", "notes"], summary_rows)
    (OUT_DIR / "round7_tree_handoff.md").write_text(
        render_handoff(summary_rows, target_presence_rows, target_neighborhood_rows, manual_review_rows, flagged_neighborhood_rows),
        encoding="utf-8",
    )
    manifest_rows = build_manifest()
    write_tsv(
        "results/06_tree_len255/round7_artifact_manifest.tsv",
        [
            "relative_path",
            "artifact_role",
            "file_size_bytes",
            "md5",
            "tracked_by_git",
            "ignored_by_git",
            "required_for_next_stage",
            "source_inputs",
            "recreate_command_or_source",
            "notes",
        ],
        manifest_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
