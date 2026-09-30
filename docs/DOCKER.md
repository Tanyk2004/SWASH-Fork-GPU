# Running everything in a container (no sudo needed)

The images hold only the toolchain; the repository is mounted at `/work`, so
builds and outputs land in your checkout as usual (`build-gnu/`, `build-nvhpc/`,
`tests/work/`, `tests/golden/`).

## Glossary

* **NVIDIA HPC SDK**: NVIDIA's compiler suite; `nvfortran` in it compiles the OpenACC directives for the GPU.
* **NVIDIA Container Toolkit**: the host-side piece that lets a container see the GPU (`--gpus all`). Installed once by an administrator; it is not needed for the CPU image.
* **Rootless Docker / Podman**: container engines that run entirely as your user, without sudo.

## Prerequisites on the workstation

One of these, none of which needs sudo at run time:

* Docker with your user in the `docker` group, or rootless Docker (`dockerd-rootless-setuptool.sh install`, needs a one-time admin step on some distributions);
* Podman (rootless by default). GPU access under Podman needs the CDI spec from the Container Toolkit (`nvidia-ctk cdi generate`, an admin step) and works with `--device nvidia.com/gpu=all`.

For the GPU image the host needs an NVIDIA driver that supports the CUDA version in the image tag: 12.9 needs driver 525 or newer, 13.x needs 580 or newer. Pick the tag accordingly (see below).

**RTX 50xx (Blackwell, compute capability 120):** supported by HPC SDK 25.3 and newer, so the default tag works; the host driver must be 570 or newer (the 5090 ships with such drivers). Build with `GPU_CC=120`.

## Build the images

```bash
docker build -f docker/Dockerfile --target cpu -t swash:cpu .
docker build -f docker/Dockerfile --target gpu -t swash:gpu .          # default HPC SDK 25.7 / CUDA 12.9
docker build -f docker/Dockerfile --target gpu -t swash:gpu --build-arg NVHPC_TAG=26.1-devel-cuda13.1-ubuntu24.04 .
```

The GPU base image is about 10 GB; the pull is the slow part. Available tags are listed on the [NGC catalog](https://catalog.ngc.nvidia.com/orgs/nvidia/containers/nvhpc/tags) (for example `25.7-devel-cuda12.9-ubuntu24.04`, `25.9-devel-cuda13.0-ubuntu24.04`, `26.1-devel-cuda13.1-ubuntu24.04`).

## Run

```bash
docker/run.sh cpu docker/build-all.sh cpu        # gfortran build, goldens for the 3 cases, branch test
docker/run.sh gpu docker/build-all.sh gpu 120    # OpenACC build for an RTX 5090 (compute capability 120), GPU-tier regression
docker/run.sh gpu tests/run_case.sh build-nvhpc marina_lay2 --tier gpu
docker/run.sh cpu                                # interactive shell in the container
```

`docker/run.sh` picks docker or podman, adds the GPU flags for the `gpu`
target, mounts the repository at `/work` and keeps created files owned by
you. `docker compose -f docker/docker-compose.yml run --rm gpu ...` is the
equivalent with Compose.

Any script in the repository can be run the same way, for example the dataset
farm:

```bash
docker/run.sh gpu python3 tools/swash_farm.py --exe build-nvhpc/bin/swash.exe \
    --template tests/cases/marina_lay2/INPUT --params params.csv \
    --gen "python3 tests/cases/gen_inputs.py marina_lay2" --jobs 4 --gpus 0 --out runs/set01
```

## If the machine has no container engine at all

Apptainer (formerly Singularity) is the usual answer on shared workstations
and clusters and runs without root:

```bash
apptainer build swash-gpu.sif docker://swash:gpu          # or docker-daemon://swash:gpu after a local docker build
apptainer exec --nv --bind "$PWD":/work --pwd /work swash-gpu.sif docker/build-all.sh gpu 89
```

`--nv` exposes the host GPU driver; no Container Toolkit is needed.

## Notes

* MPI runs inside the container: `docker/run.sh cpu tests/run_case.sh build-gnu marina_lay2 --np 4` after configuring with `-DMPI=ON`. The CPU image ships OpenMPI; the GPU image uses the HPC SDK's own OpenMPI.
* The container runs as your user id (docker) or maps root to you (podman), so `git` inside the container sees your checkout with the right ownership.
* Generated `.f90` files and build directories are excluded from the image by `.dockerignore`; they only ever live in the mounted checkout.
