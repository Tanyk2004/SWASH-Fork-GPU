# Building and validating the GPU / engine fork

Everything below runs on your CUDA workstation. The CPU reference build needs
only gfortran; the GPU build needs the NVIDIA HPC SDK (`nvfortran`).

## Glossary

* **OpenACC**: `!$acc` directives that offload loops to the GPU when compiled with `nvfortran`; plain comments for gfortran.
* **TIMG**: SWASH's built-in section timers, printed at the end of the `PRINT` file.
* **Golden**: the reference outputs of a case produced by the trusted CPU build.

## 1. CPU reference build (gfortran)

```bash
mkdir build-gnu && cd build-gnu
cmake .. -G Ninja -DTIMG=ON -DENGINE=ON
cmake --build .
cd ..
```

Then create the goldens and the profiling baseline:

```bash
tests/run_case.sh build-gnu marina_lay2  --golden     # layered non-hydrostatic, explicit (GPU target path)
tests/run_case.sh build-gnu basin_da     --golden     # depth-averaged, SIP solver
tests/run_case.sh build-gnu marina_float --golden     # fixed pontoon, implicit scheme
```

Each run prints the timing table (`tools/profile_timg.py`); copy the numbers
into `docs/PERF.md`. The goldens land in `tests/golden/<case>/`.

To confirm that the driver refactor (`src/SwashRun.ftn90`) is behaviour
preserving, build the pristine upstream tree in a second directory and compare
its outputs with `tests/regress.py --tier cpu --atol 0 --rtol 0` (bitwise).

## 2. Engine checks

```bash
python3 tests/test_branch.py build-gnu/lib/libswash_engine.so marina_lay2
```

must print `PASS` (save / restore / re-step is bitwise identical). A quick
interactive smoke test:

```python
import sys; sys.path.insert(0, "engine/python")
from swash_engine import Swash, hull_pressure
sw = Swash("build-gnu/lib/libswash_engine.so")
sw.init("tests/work/marina_lay2-build-gnu")       # a directory prepared by run_case.sh
sw.step(100)
eta = sw.field("WATL", masked=True)
sw.set_pressure(hull_pressure(sw, x=60, y=75, heading=0.0, length=30, beam=8, draft=1.5))
sw.step(500)                                       # a Kelvin-type wake develops behind the patch if you move it each step
sw.finish()
```

## 3. GPU build (nvfortran, OpenACC)

```bash
mkdir build-nvhpc && cd build-nvhpc
cmake .. -G Ninja -DCMAKE_Fortran_COMPILER=nvfortran -DOPENACC=ON -DGPU_CC=89 -DTIMG=ON -DENGINE=ON
cmake --build . 2>&1 | tee build.log
cd ..
```

`GPU_CC` is the compute capability: 80 (A100), 86 (RTX 30xx), 89 (RTX 40xx),
90 (H100). Check `build.log` for the `-Minfo=accel` messages of
`SwashExpLay2DHflow.f90`: every ported loop must report "Generating NVIDIA GPU
code" and the collapsed loops "gang vector"; a "loop not parallelizable" line
on one of them is a bug to fix before running.

Validate against the goldens with the physics tier:

```bash
tests/run_case.sh build-nvhpc marina_lay2 --tier gpu
NV_ACC_NOTIFY=3 tests/run_case.sh build-nvhpc marina_lay2 --tier gpu 2>&1 | grep -c 'upload\|download'
```

The second command counts data transfers; during the incremental port each
ported region uploads its inputs and downloads its outputs on purpose (see
`src/SwashAccData.ftn90`), so the count is large. It must drop to a handful per
step once the residency is complete (Phase 1 done).

CPU fallback of the same binary: run with `ACC_DEVICE_TYPE=host` (kernels on
the CPU) to separate directive bugs from GPU-specific ones.

## 4. Where the port stands (first milestone)

Ported regions of `src/SwashExpLay2DHflow.ftn90`:

| region | lines (approx.) | pattern |
|---|---|---|
| maximum CFL number | `[ACC region: CFL]` | 3-level collapse, `reduction(max:)` |
| global continuity (mass fluxes, water level) | `[ACC region: continuity]` | 2/3-level collapse, sequential layer loops |
| vertical terms + tridiagonal solve of u-momentum | `[ACC region: vertical u-momentum]` | one thread per column, `fluxlim_dev` |

Next regions in order: v-momentum vertical terms (mirror of u), w-momentum,
Poisson matrix build, pressure update and velocity correction, the explicit
horizontal terms of u and v (split at the halo exchanges), sponge, then the
routines around the kernel (`SwashUpdateDepths`, `SwashLayerIntfaces`,
`SwashDryWet`, `SwashBreakPoint`). The pressure solve itself gets a GPU
BiCGSTAB with a column preconditioner (`src/SwashSolversAcc.ftn90`, Phase 1b).

## 5. Parallel dataset generation

```bash
python3 tools/swash_farm.py --exe build-gnu/bin/swash.exe \
    --template tests/cases/marina_lay2/INPUT --params params.csv \
    --gen "python3 tests/cases/gen_inputs.py marina_lay2" --jobs 4 --gpus 0,1 --out runs/set01
```

`params.csv` columns are substituted into `{placeholders}` in the template
(make a copy of the deck with e.g. `SPECTRUM {hs} {tp} {dir} 20.` and
`SEED {seed}`). With GPU builds, start `nvidia-cuda-mps-control -d` first so
several processes share a GPU efficiently.

## 6. Notes and known limits

* `switch.pl` (the `!MPI`/`!TIMG`/`!ACC` line preprocessor) only regenerates a
  `.f90` when its `.ftn90` is newer, so after changing build options run
  `cmake -P clobber.cmake` or delete `src/*.f90` before reconfiguring.
* The engine is single-process (one SWASH instance per process, single MPI rank).
* Bathymetry changes through the engine re-derive the velocity-point depths;
  a cell that dries keeps its last velocity until the wet/dry masks update on
  the next step (same as SWASH's own behaviour after a hot start).
* Boundary values read from `SERIES` files are not rewound by `restore`; use
  generated spectra or regular waves for branch-and-predict runs.
