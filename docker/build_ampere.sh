#!/usr/bin/env bash
# Ampere (sm_80) — default public image (A100 / many cloud GPUs).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
exec docker build -f docker/Dockerfile -t 4dwcm:ampere \
  --build-arg CUDA_ARCHITECTURES=80 \
  --build-arg GPU_ARCH=sm_80 \
  --build-arg KOKKOS_ARCH_FLAGS="-DKokkos_ARCH_AMPERE80=yes" \
  --build-arg BTREE_NVCC_ARCH=sm_80 \
  "$@" \
  .
