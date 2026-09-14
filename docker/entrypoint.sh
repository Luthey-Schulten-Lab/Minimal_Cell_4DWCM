#!/usr/bin/env bash
# Activate the Lattice Microbes conda env and set DNA software defaults.
set -euo pipefail
source /opt/conda/etc/profile.d/conda.sh
conda activate lm_2.5_dev
export PATH="/Software/sc_chain_generation/src:/Software/btree_chromo/build/apps:${PATH}"
export PYTHONPATH="/Software/Lattice_Microbes/src/pylm:/Software/Lattice_Microbes/src/jlm:${PYTHONPATH:-}"
cd /src/4d
exec "$@"
