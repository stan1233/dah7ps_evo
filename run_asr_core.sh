#!/bin/bash
# Run IQ-TREE ASR on the pruned ingroup tree
# Output will be generated in results/04_phylogeny_asr/ASR_core.*

cd /home/tynan/0218
source /home/tynan/miniforge3/etc/profile.d/conda.sh
conda activate dah7ps_v4

iqtree -s results/03_msa_core/core_asr.afa \
  -te results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
  -m Q.PFAM+F+R10 \
  -asr \
  -T 20 \
  --prefix results/04_phylogeny_asr/ASR_core
