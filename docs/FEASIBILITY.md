# Feasibility of a GPU-accelerated, real-time, coupled SWASH

Assessment of the SWASH 12.01 source in this repository against the goals of the
fork: GPU acceleration without changing the physics, parallel instances for
synthetic datasets, real-time environment editing, two-way coupling with
Gazebo and Isaac Sim (vessel wakes), a future predictor, and a 3D viewer.

## Glossary

| Term | Meaning |
|---|---|
| SWASH | Simulating WAves till SHore, the non-hydrostatic wave-flow model (~180k lines Fortran 90) |
| 2DH | two-dimensional horizontal grid, with optional vertical layers |
| non-hydrostatic pressure `q` | the pressure term that gives correct wave dispersion; solved every step from a Poisson-type system |
| Keller-box | the default vertical discretisation of `q` |
| BiCGSTAB | Bi-Conjugate Gradient Stabilised, the iterative solver of the pressure system |
| ILU / ILUD | Incomplete LU factorisation preconditioners; their triangular sweeps are sequential |
| SIP | Stone's Implicit Procedure, the solver of 5-point (depth-averaged) systems |
| CFL | Courant-Friedrichs-Lewy number; the explicit scheme keeps it in [0.4, 0.8] by halving/doubling `dt` |
| OpenACC | `!$acc` compiler directives for GPU offload, compiled with `nvfortran` (NVIDIA HPC SDK); comments to other compilers |
| MPI | Message Passing Interface, SWASH's existing multi-process parallelism |
| `kgrpnt(m,n)` | index map from 2D cell (m,n) to the compressed 1D point index `nm`; `nm = 1` is a dummy slot for inactive cells |
| `patm` | atmospheric-pressure field; its gradient enters the momentum equations; also usable as a vessel's pressure footprint |
| Kelvin wake | the V-shaped wave pattern behind a moving vessel; produced by a translating pressure patch |
| FSI | fluid-structure interaction |

## Verdict per use case

| Use case | Verdict | Evidence in the source |
|---|---|---|
| Many parallel instances for datasets | **Feasible today, one instance per process** | All state is module-global `save, allocatable` (`src/SwashFlowdata.ftn90`, `src/SwashModule1.ftn90`, `src/SwashModule2.ftn90`); fixed unit numbers 3/4 and fixed file names `INPUT`, `PRINT`, `Errfile`, `norm_end` in the working directory (`src/ocpids.ftn:122-304`); `SAVE` entry counters in 146 files. Several instances in one process are impossible; N processes in N directories work (`tools/swash_farm.py`). |
| Real-time environment editing | **Feasible with the engine layer** | Current, friction, wind, water level and pressure fields are re-read every step (`src/SwashUpdateData.ftn90:2489-2518`); the pressure gradient enters every solver (e.g. `src/SwashExpLay2DHflow.ftn90:1145`). Bathymetry `dps` is fixed within a COMPUTE block; the engine re-derives the velocity-point depths after a change (`SwashFlowDP`). The draft (vessel) field is not updated at run time, so vessels go through the pressure field. |
| Two-way Gazebo / Isaac Sim coupling | **Feasible as pose-in / surface-and-forces-out; hull-in-water FSI later** | A moving pressure patch (`swash_set_pressure_field`) creates wakes without touching the discretisation. Water forces on a hull follow from integrating hydrostatic + non-hydrostatic pressure under the footprint (the same integral as `src/SwashHydroLoads.ftn90`). SWASH's own floating bodies (`ifloat=2`, 6 degrees of freedom, `src/SwashMotionRigidBod.ftn90`) are small-motion around a fixed footprint and require the implicit scheme (`src/SwashCheckPrep.ftn90:1113-1150`; footprint mapped once at `:3148`), so a translating hull is out of first scope. |
| Future predictor (branch and run ahead) | **Feasible with the in-memory checkpoint** | The hot file is ASCII, written only between COMPUTE blocks and omits `q`, `w`, turbulence and body state (`src/SwashBackup.ftn90:101-363`, single call site `src/SwashReadInput.ftn90:2817`). All state is a fixed set of arrays allocated once per COMPUTE (`SwashCheckPrep`, 131 allocates), so `engine/swash_state.ftn90` (generated) copies every array of `SwashFlowdata` and `SwashRigBoddata` plus the time scalars. Running ahead faster than real time depends on the GPU speed-up (see `docs/PERF.md`). |
| 3D viewer | **Feasible** | VTK output exists but writes a flat z = 0 surface (`src/SwashVTKWriteData.ftn90:219-229`); the engine gives a live path to the surface and layer data. |

## GPU port: feasible, one hard part

* Nearly every per-step operation is a pointwise or stencil loop over `kgrpnt`
  with a short vertical loop per column. The vertical implicit solves are
  already one independent tridiagonal system per column
  (`src/SwashExpLay2DHflow.ftn90:1947-1970`).
* Precision is 32-bit `REAL` (no `-r8` flag); arrays are allocated once per
  COMPUTE; no allocation inside kernels.
* Per-step call sequence (structured grid): `SwashUpdateData` (boundary values,
  forcing), `SwashComputStruc` (momentum, pressure Poisson, continuity, then
  depths, layers, wet/dry, breaking, transport, turbulence), `SwashOutput`,
  adaptive `dt` (`src/SwashMain.ftn90`, now `src/SwashRun.ftn90`).
* **The hard part is the pressure solver's preconditioner.** BiCGSTAB itself
  (`src/SwashSolvers.ftn90:3312`) is GPU-friendly (stencil matrix-vector
  products, dot products), but the ILU preconditioners (`ilu`, `iludr`,
  `iluds`, applied in `ivl`/`ivu`, `src/SwashSolvers.ftn90:610-935`) are
  lexicographic recurrences rebuilt every step. They must be replaced on the
  GPU by a parallel preconditioner (block-Jacobi per water column first). Under
  MPI the ILU is already block-Jacobi per subdomain, so the converged solution
  already depends on the partitioning only to solver tolerance; the same holds
  for the replacement.
* Index order: points are numbered y-fastest (`src/SwashInitCompGrid.ftn90:274-282`)
  while loops run with x innermost. The port keeps the numbering and collapses
  the loops so consecutive threads get consecutive points.
* `fluxlim` (`src/SwashServices.ftn90:307`) reads module state and can call
  the error logger; a pure device copy `fluxlim_dev` is used inside kernels.

## Prior art

* FUNWAVE-TVD GPU (CUDA Fortran): 4-7x over a 36-core node.
* Celeris (Boussinesq, GPU): faster than real time, interactive editing.
* XBeach on GPU (NIWA): desktop GPUs competitive with CPU clusters.
* asv_wave_sim (Gazebo): wave visual and buoyancy plugins driven by one wavefield topic; the adapter design follows it.

## What this fork adds (first milestone)

* Build options `OPENACC`, `GPU_CC`, `TIMG`, `ENGINE` and `nvfortran` support (`cmake/SwashOptions.cmake`, `src/CMakeLists.txt`).
* Reference cases and regression harness (`tests/`), per-section timers (`src/SwashTimers.ftn90`, `tools/profile_timg.py`).
* OpenACC data residency (`src/SwashAccData.ftn90`) and the first ported regions of `src/SwashExpLay2DHflow.ftn90` (CFL reduction, continuity, vertical u-momentum solve).
* Re-entrant driver (`src/SwashRun.ftn90`) and the engine API (`engine/`), including in-memory checkpoints and the external pressure / mass-source / bathymetry hooks (`src/SwashExtData.ftn90`).
