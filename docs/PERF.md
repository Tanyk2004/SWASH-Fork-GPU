# Performance log

Fill this in from runs on the CUDA workstation (see `docs/BUILD_GPU.md`).
Each row: case, build, hardware, wall time of the compute phase, step time,
and the split from `tools/profile_timg.py`.

## Baseline (CPU, gfortran, 1 core)

| case | grid | steps | compute wall [s] | ms/step | flow solver % | ILU build % | BiCGSTAB % | notes |
|---|---|---|---|---|---|---|---|---|
| marina_lay2 | 200x150x2 | | | | | | | |
| marina_float | 200x150x2 | | | | | | | implicit, ILU (icond=3) |
| basin_da | 300x300x1 | | | | | | | SIP |

## GPU (nvfortran, OpenACC)

| case | GPU | ms/step | speed-up vs 1 core | real-time factor | iterations/step (CPU -> GPU) | notes |
|---|---|---|---|---|---|---|

## Expectations (to be confirmed)

| Grid | CPU 1 core / step | CPU 16 MPI ranks | GPU (A100 / RTX 4090) | real-time factor on GPU |
|---|---|---|---|---|
| 200x150x2 | 25-40 ms | 4-6 ms | 2-3 ms | 20-30x |
| 500x500x3 | 0.4-0.6 s | 50-80 ms | 12-25 ms | 3-6x |
| 1000x1000x3 | 2-3 s | 0.3-0.4 s | 50-90 ms | ~1x |

Basis: FUNWAVE-GPU reports 4-7x over a 36-core node, Celeris runs Boussinesq
faster than real time; a layered non-hydrostatic step is 3-4x costlier per
cell because of the pressure solve. The estimates assume the pressure solve on
the GPU with a column preconditioner at about twice the ILU iteration count.
