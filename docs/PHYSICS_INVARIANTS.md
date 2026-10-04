# Physics invariants of the GPU port

The purpose of the fork is speed and coupling, not a different model. This
file states what the port is allowed to change and what it is not, and how
that is checked.

## Allowed to change

* **Loop nesting and order** of point loops (each grid point is independent
  inside a kernel).
* **Order of reductions**: the maximum CFL number, the dot products and norms
  inside BiCGSTAB, mass/energy diagnostics. This changes round-off only.
* **The preconditioner** of the pressure solver (ILU family -> column
  block-Jacobi or diagonal on the GPU). It changes the number of iterations,
  not the equations; the converged pressure is the same to the solver
  tolerance. Under MPI the ILU is already subdomain-local, so this class of
  change is already present in the reference code.
* **Placement of halo exchanges** and host/device copies.
* **The implementation form** of `fluxlim` (`fluxlim_dev` is the same formulas
  without module state, tracing or error I/O).
* **Writes to the dummy point** `nm = 1` may be skipped in kernels; every
  consumer either resets it or masks it (`uwetp(1) = 0`, `s1(1) = 0`).

## Not allowed to change

* Any discretisation coefficient, stencil, limiter formula, or the order of
  the schemes (MacCormack corrections, Keller-box, Hancock/RK3 continuity).
* The adaptive time-step rule (`SwashRunStep`: halve above `CFLHIG`, double
  below `CFLLOW` after 20 constant steps).
* Solver tolerance defaults (`RHSACCUR`, `MAXITER`), the default
  preconditioner of the CPU path (`icond = 2`), wet/dry thresholds, breaking
  criteria.
* Boundary condition generation (spectra, Fourier components, random phases).
* Output quantity definitions and file formats.
* The external forcing hooks add the **same terms through the same code** as
  the file-driven inputs: `patm` (pressure gradient), `srcm` (mass source),
  `dps` (bathymetry). They never add new terms.

## How it is checked

| Check | Tool | Tolerance |
|---|---|---|
| Refactor of the driver (`SwashRun`) | `tests/run_case.sh` on the pristine and the forked build | bitwise identical outputs |
| CPU compilers / flags | `tests/regress.py --tier cpu` | `atol 1e-4 m`, `rtol 1e-4`; gauge Hm0 within 0.5 %, Tp within 1 % |
| GPU vs CPU, same compiler (nvfortran CPU build as baseline; cross-compiler runs are different random wave realisations because of the intrinsic `random_number`, so only their gauge statistics are comparable) | `tests/regress.py --tier gpu` | `atol 1e-3 m`, `rtol 1e-3`; Hm0 within 1 %, Tp within 2 % |
| GPU preconditioner | same case with `RHSACCUR 1e-6` on both | water level within `1e-5 m` |
| Volume conservation | `sum(gsqs * (s1 + dps))` drift per run, both builds | relative drift below `1e-6` |
| Checkpoint completeness | `tests/test_branch.py` | bitwise identical after save/restore |
| Data movement | `NV_ACC_NOTIFY=3` on a GPU run | no implicit transfers inside the time loop (once residency is complete) |
