#!/usr/bin/env bash
# Run one reference case with a given build and (optionally) compare against goldens.
#
#   tests/run_case.sh <build-dir> <case> [--golden] [--compare <golden-dir>] [--tier cpu|gpu] [--np N]
#
#   <build-dir>   directory holding bin/swash.exe (e.g. build-gnu, build-nvhpc)
#   <case>        name under tests/cases/
#   --golden      copy the outputs into tests/golden/<case>/ (creates the reference)
#   --compare D   run tests/regress.py against golden directory D (default tests/golden/<case>)
#   --tier        tolerance tier for the comparison (cpu = CPU-vs-CPU, gpu = physics tier)
#   --np N        run under mpirun -np N (build must have MPI=ON)
#
# Outputs land in tests/work/<case>-<build-name>/ together with PRINT and the timing table.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/.." && pwd)"

build="$1"; shift
case_name="$1"; shift
golden=0; compare=""; tier="cpu"; np=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --golden) golden=1;;
    --compare) compare="$2"; shift;;
    --tier) tier="$2"; shift;;
    --np) np="$2"; shift;;
    *) echo "unknown option $1" >&2; exit 2;;
  esac
  shift
done

exe="$(cd "$build" && pwd)/bin/swash.exe"
[[ -x "$exe" ]] || { echo "no executable at $exe" >&2; exit 1; }
[[ -f "$here/cases/$case_name/INPUT" ]] || { echo "no case $case_name" >&2; exit 1; }

work="$here/work/${case_name}-$(basename "$build")"
rm -rf "$work"; mkdir -p "$work"
cp "$here/cases/$case_name/INPUT" "$work/INPUT"
python3 "$here/cases/gen_inputs.py" "$case_name" "$work"

echo "== running $case_name with $exe in $work"
start=$(date +%s.%N)
if [[ "$np" -gt 0 ]]; then
  ( cd "$work" && mpirun -np "$np" "$exe" )
else
  ( cd "$work" && "$exe" )
fi
end=$(date +%s.%N)
echo "== wall time: $(python3 -c "print(f'{$end-$start:.1f} s')")"

if [[ -f "$work/norm_end" ]]; then echo "== normal end"; else echo "== NO norm_end: check $work/PRINT and Errfile" >&2; exit 1; fi
if [[ -f "$work/Errfile" ]]; then echo "== Errfile present:"; cat "$work/Errfile"; fi
grep -i 'warning\|error' "$work/PRINT" | head -20 || true

# timing table (only present for TIMG builds)
if grep -q '# Details on timings' "$work/PRINT"; then
  python3 "$root/tools/profile_timg.py" "$work/PRINT" | tee "$work/timing.txt"
fi

if [[ "$golden" -eq 1 ]]; then
  dest="$here/golden/$case_name"
  mkdir -p "$dest"
  cp "$work"/*.mat "$work"/*.tbl "$dest"/ 2>/dev/null || true
  cp "$work/PRINT" "$dest/PRINT"
  echo "== goldens written to $dest"
fi

if [[ -n "$compare" || "$golden" -eq 0 ]]; then
  ref="${compare:-$here/golden/$case_name}"
  if [[ -d "$ref" ]]; then
    python3 "$here/regress.py" "$ref" "$work" --tier "$tier"
  else
    echo "== no golden directory $ref; run with --golden on the reference CPU build first"
  fi
fi
