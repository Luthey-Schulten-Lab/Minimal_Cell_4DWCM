#!/usr/bin/env bash
# Short smoke run inside a built 4dwcm image (requires NVIDIA Container Toolkit).
set -euo pipefail
IMAGE="${1:-4dwcm:ampere}"
OUT="${2:-docker_smoke}"
BIO_T="${3:-60}"

docker run --rm --gpus all \
  -e DNA_GPU_ID=0 \
  -v "${PWD}/Data:/src/4d/Data" \
  "${IMAGE}" \
  python Whole_Cell_Minimal_Cell.py \
    -od "${OUT}" \
    -t "${BIO_T}" \
    -cd 0 \
    -drs 1 \
    -dsd /Software/
