#!/usr/bin/env bash
# Inside the container: configure, build, and validate.
#
#   docker/build-all.sh cpu            gfortran build (build-gnu), goldens for the three cases, branch test
#   docker/build-all.sh gpu [GPU_CC]   nvfortran/OpenACC build (build-nvhpc) and the GPU-tier regression
#                                      GPU_CC: 80 A100, 86 RTX 30xx, 89 RTX 40xx, 90 H100 (default: native)
set -euo pipefail
target="${1:-cpu}"; cc="${2:-native}"
cd "$(dirname "$0")/.."
case "$target" in
  cpu)
    rm -rf build-gnu src/*.f90 src/*.f engine/*.f90
    cmake -S . -B build-gnu -G Ninja -DTIMG=ON -DENGINE=ON
    cmake --build build-gnu
    for c in marina_lay2 basin_da marina_float; do tests/run_case.sh build-gnu "$c" --golden; done
    python3 tests/test_branch.py build-gnu/lib/libswash_engine.so marina_lay2
    ;;
  gpu)
    command -v nvfortran >/dev/null || { echo "nvfortran not found: use the gpu image" >&2; exit 1; }
    rm -rf build-nvhpc src/*.f90 src/*.f engine/*.f90
    cmake -S . -B build-nvhpc -G Ninja -DCMAKE_Fortran_COMPILER=nvfortran -DOPENACC=ON -DGPU_CC="$cc" -DTIMG=ON -DENGINE=ON
    cmake --build build-nvhpc 2>&1 | tee build-nvhpc/build.log
    echo "== OpenACC report for the ported kernel:"
    grep -A2 'SwashExpLay2DHflow\|swashexplay2dhflow' build-nvhpc/build.log | grep -i 'gpu code\|not parallel\|gang' | head -40 || true
    [[ -d tests/golden/marina_lay2 ]] || { echo "no goldens yet: run docker/build-all.sh cpu first"; exit 1; }
    tests/run_case.sh build-nvhpc marina_lay2 --tier gpu
    ;;
  *) echo "usage: $0 cpu|gpu [GPU_CC]" >&2; exit 2;;
esac
