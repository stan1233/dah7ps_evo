#!/usr/bin/env bash
set -euo pipefail

# Formal S1-only KDOPS11 rerun after excluding the KDOPS_O66496 outlier.
# This intentionally does not run S2. The output must remain distinct from the
# legacy 12-KDOPS S1/S2 artifacts.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

IQTREE_BIN="${IQTREE_BIN:-/home/luogu/.local/bin/iqtree3}"
PYTHON_BIN="${PYTHON_BIN:-/home/luogu/miniforge3/envs/dah7ps_v4/bin/python}"
THREADS="${THREADS:-20}"

ALIGNMENT="results/04_phylogeny_asr/core_with_outgroup_KDOPS11_noO66496.afa"
PREFIX="results/04_phylogeny_asr/CoreTree_rooted_MFP_KDOPS11_noO66496"
TREE="${PREFIX}.treefile"
INGROUP_TREE="${PREFIX}_ingroup.treefile"
COMPARE_TSV="results/04_phylogeny_asr/S1_KDOPS11_noO66496_vs_legacy12_ingroup_comparison.tsv"
COMPARE_MD="results/04_phylogeny_asr/S1_KDOPS11_noO66496_vs_legacy12_ingroup_comparison.md"

if [[ ! -x "$IQTREE_BIN" ]]; then
  echo "ERROR: IQ-TREE executable not found or not executable: $IQTREE_BIN" >&2
  exit 1
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "ERROR: Python executable not found or not executable: $PYTHON_BIN" >&2
  exit 1
fi

if [[ ! -f "$ALIGNMENT" ]]; then
  echo "ERROR: KDOPS11 alignment not found: $ALIGNMENT" >&2
  exit 1
fi

for output in "$TREE" "$INGROUP_TREE" "$COMPARE_TSV" "$COMPARE_MD"; do
  if [[ -e "$output" ]]; then
    echo "ERROR: output already exists; refusing to overwrite: $output" >&2
    exit 1
  fi
done

"$IQTREE_BIN" \
  -s "$ALIGNMENT" \
  -m MFP \
  -B 1000 \
  -T "$THREADS" \
  -o KDOPS_P0A715,KDOPS_Q9ZFK4 \
  --prefix "$PREFIX"

"$PYTHON_BIN" scripts/prune_tree.py \
  --input "$TREE" \
  --remove_prefix KDOPS_ \
  --output "$INGROUP_TREE" \
  --assert_rooted

"$PYTHON_BIN" scripts/assert_tip_match.py \
  --tree "$INGROUP_TREE" \
  --msa results/03_msa_core/core_asr.afa \
  --assert_identical

"$PYTHON_BIN" scripts/compare_rooted_tree_pairs.py \
  --left results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
  --right "$INGROUP_TREE" \
  --left-label S1_legacy12_pruned_ingroup \
  --right-label S1_KDOPS11_noO66496_ingroup \
  --output-tsv "$COMPARE_TSV" \
  --output-md "$COMPARE_MD" \
  --membership-fasta Ia=results/02_qc/nr80_Ia.fasta \
  --membership-fasta Ib=results/02_qc/nr80_Ib.fasta \
  --membership-fasta II=results/02_qc/nr80_II.fasta
