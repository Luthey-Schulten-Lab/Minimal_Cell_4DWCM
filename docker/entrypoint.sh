#!/usr/bin/env bash
# Activate the Lattice Microbes conda env and set DNA software defaults.
set -euo pipefail
# Conda's activate-binutils_linux-64.sh evaluates $HOST, which is unset, so
# nounset aborts activation. Lift -u across the activate call only; -e stays
# on, so a genuine activation failure still fails the container.
set +u
source /opt/conda/etc/profile.d/conda.sh
conda activate lm_2.5_dev
set -u
export PATH="/Software/sc_chain_generation/src:/Software/btree_chromo/build/apps:${PATH}"
# pylm and jLM are installed into the conda env by `make install`; adding the
# source trees here would shadow the built jLM egg (see Dockerfile note).
export PYTHONPATH="${PYTHONPATH:-}"
cd /src/4d
exec "$@"
