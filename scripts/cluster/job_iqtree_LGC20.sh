#!/bin/bash

# --- JH Unischeduler Options ---
#JSUB -q normal
#JSUB -n 60
#JSUB -R "span[hosts=1]"
#JSUB -M 128000
#JSUB -J "iqtree_core_LGC20"
#JSUB -e log/error_LGC20.%J
#JSUB -o log/output_LGC20.%J

# ============================================================
# IQ-TREE: LG+C20+F+G rooted tree with KDOPS outgroup
#
# Purpose : Phase 4.1 S2 scenario — site-heterogeneous model
#           for root robustness comparison with MFP (S1).
# Input   : core_with_outgroup.afa (9,405 seqs × 472 cols)
# Output  : CoreTree_rooted_LGC20.*
#
# Notes   :
#   - LG+C20+F+G requires ~55 GB RAM; request 128 GB for safety.
#   - Local run was OOM-killed; must run on cluster.
#   - Use checkpoint resume (-ckp) if restarting from a crash.
#   - Outgroup: 12 KDOPS sequences.
# ============================================================

echo "=== Job start ==="
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "Cores: $LSB_DJOB_NUMPROC"
echo ""

# --- Paths ---
WORKDIR=/gpfshddpool/home/luoguangyuan/ssd/0309_iqtree
cd "$WORKDIR" || { echo "ERROR: cannot cd to $WORKDIR"; exit 1; }

mkdir -p log

# --- Run IQ-TREE ---
~/mambaforge/envs/iqtree3/bin/iqtree \
    -s core_with_outgroup.afa \
    -m LG+C20+F+G \
    -B 1000 \
    -T 60 \
    -o KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496 \
    --prefix CoreTree_rooted_LGC20

EXIT_CODE=$?

echo ""
echo "=== Job end ==="
echo "Date: $(date)"
echo "Exit code: $EXIT_CODE"

exit $EXIT_CODE
