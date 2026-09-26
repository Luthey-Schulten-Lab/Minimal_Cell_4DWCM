#!/bin/bash
# Build py/bdry/bdry_gen: btree_chromo's own boundary_surface.cpp behind a tiny CLI (the membrane shell's vertex order comes from
# libstdc++'s unordered_set iteration, so it must be btree_chromo's code built with the same compiler as btree_chromo itself).
#   usage: build.sh <btree_chromo source dir>        env: WCM_IMAGE (build inside that image, recommended), CXX (default g++)
set -eu; D=$(dirname "$(readlink -f "$0")"); SRC=$(realpath "${1:?btree_chromo source dir}")
CMD="\${CXX:-g++} -O2 -std=c++17 -I$SRC/include $D/bdry_gen.cpp $SRC/src/boundary_surface.cpp $SRC/src/vec_quat_manipulator.cpp -lm -o $D/bdry_gen && \${CXX:-g++} --version | head -1"
if [ -n "${WCM_IMAGE:-}" ]; then
  docker run --rm --user "$(id -u):$(id -g)" -v "$SRC:$SRC:ro" -v "$D:$D" --entrypoint bash "$WCM_IMAGE" -c "export PATH=/usr/local/bin:\$PATH; $CMD"
else bash -c "$CMD"; fi
echo "built $D/bdry_gen"
