#!/usr/bin/env bash
# Blackwell (sm_100) — NVIDIA B200 production box.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
exec docker build -f docker/Dockerfile -t 4dwcm:blackwell \
  --build-arg CUDA_ARCHITECTURES=100 \
  --build-arg GPU_ARCH=sm_100 \
  --build-arg KOKKOS_ARCH_FLAGS="-DKokkos_ARCH_BLACKWELL100=yes" \
  --build-arg BTREE_NVCC_ARCH=sm_100 \
  "$@" \
  .
