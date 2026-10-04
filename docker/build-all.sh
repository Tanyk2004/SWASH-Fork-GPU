#!/usr/bin/env bash
# Inside the container: configure, build, and validate.
#
#   docker/build-all.sh cpu            gfortran build (build-gnu), goldens for the three cases, branch test
#   docker/build-all.sh gpu [GPU_CC]   nvfortran CPU baseline (build-nvhpc-cpu), nvfortran/OpenACC build (build-nvhpc),
#                                      GPU-vs-same-compiler-CPU regression
#                                      GPU_CC: 80 A100, 86 RTX 30xx, 89 RTX 40xx, 90 H100, 120 RTX 50xx (Blackwell; default: native)
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
    # (a) same-compiler CPU baseline: nvfortran without OpenACC.  SWASH draws its random wave
    #     phases with the Fortran intrinsic random_number, which differs between compilers, so a
    #     GPU run can only be compared point by point with the same compiler's CPU run.
    rm -rf build-nvhpc-cpu src/*.f90 src/*.f engine/*.f90
    cmake -S . -B build-nvhpc-cpu -G Ninja -DCMAKE_Fortran_COMPILER=nvfortran -DOPENACC=OFF -DTIMG=ON
    cmake --build build-nvhpc-cpu 2>&1 | tail -3
    tests/run_case.sh build-nvhpc-cpu marina_lay2 --golden --golden-dir tests/golden-nvhpc/marina_lay2
    # (b) GPU build
    rm -rf build-nvhpc src/*.f90 src/*.f engine/*.f90
    cmake -S . -B build-nvhpc -G Ninja -DCMAKE_Fortran_COMPILER=nvfortran -DOPENACC=ON -DGPU_CC="$cc" -DTIMG=ON -DENGINE=ON
    set +e
    cmake --build build-nvhpc 2>&1 | tee build-nvhpc/build.log
    rc=${PIPESTATUS[0]}
    set -e
    if [[ $rc -ne 0 ]]; then
      echo "== BUILD FAILED; first compiler errors:"
      grep -n 'NVFORTRAN-S\|NVFORTRAN-F\|Error\|error:' build-nvhpc/build.log | head -20
      exit $rc
    fi
    echo "== OpenACC report for the ported kernel (from -Minfo=accel):"
    awk 'tolower($0) ~ /^swashexplay2dhflow:/ {p=1; print; next} /^[a-z_0-9]+:$/ {p=0} p' build-nvhpc/build.log \
      | grep -i 'generating\|gang\|vector\|seq\|not parallel\|reduction\|update' | head -60 || true
    # (c) the test that matters: nvfortran GPU vs nvfortran CPU
    echo "== GPU vs same-compiler CPU baseline:"
    tests/run_case.sh build-nvhpc marina_lay2 --compare tests/golden-nvhpc/marina_lay2 --tier gpu
    # (d) informational: nvfortran CPU vs gfortran goldens (different random realisations, only the
    #     gauge statistics are comparable, and a 60 s record is short)
    if [[ -d tests/golden/marina_lay2 ]]; then
      echo "== informational: nvfortran-CPU vs gfortran goldens (different wave realisations):"
      python3 tests/regress.py tests/golden/marina_lay2 tests/work/marina_lay2-build-nvhpc-cpu --tier gpu 2>/dev/null | grep 'gauge\|== ' || true
    fi
    ;;
  *) echo "usage: $0 cpu|gpu [GPU_CC]" >&2; exit 2;;
esac
