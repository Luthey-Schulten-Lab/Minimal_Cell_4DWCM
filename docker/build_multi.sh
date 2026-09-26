#!/usr/bin/env bash
# Multi-arch — longer build / larger image; wider GPU support.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
exec docker build -f docker/Dockerfile -t 4dwcm:multi \
  --build-arg CUDA_ARCHITECTURES="70;75;80;86;89;90" \
  --build-arg GPU_ARCH=sm_80 \
  --build-arg KOKKOS_ARCH_FLAGS="-DKokkos_ARCH_VOLTA70=yes -DKokkos_ARCH_TURING75=yes -DKokkos_ARCH_AMPERE80=yes -DKokkos_ARCH_AMPERE86=yes" \
  --build-arg BTREE_NVCC_ARCH="sm_70 sm_75 sm_80 sm_86 sm_89 sm_90" \
  "$@" \
  .
