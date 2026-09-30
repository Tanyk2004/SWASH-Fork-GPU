# syntax=docker/dockerfile:1
#
# Container images for the SWASH GPU / engine fork.  Two targets:
#
#   cpu  - Ubuntu 24.04 + gfortran: reference build, regression harness, engine
#   gpu  - NVIDIA HPC SDK (nvfortran, CUDA): everything above plus the OpenACC build
#
# Build (from the repository root):
#   docker build -f docker/Dockerfile --target cpu -t swash:cpu .
#   docker build -f docker/Dockerfile --target gpu -t swash:gpu .
#
# The images contain the toolchain only; the source tree is mounted at /work
# at run time (see docker/run.sh), so edits on the host are used directly.

ARG NVHPC_TAG=25.7-devel-cuda12.9-ubuntu24.04

# ---------------------------------------------------------------- common layer
FROM ubuntu:24.04 AS cpu
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran gcc g++ cmake ninja-build make perl git \
        python3 python3-pip python3-venv \
        libopenmpi-dev openmpi-bin \
        ca-certificates curl less vim-tiny \
    && rm -rf /var/lib/apt/lists/*
# python tools used by tests/ and tools/ (numpy, scipy for regress.py; pyzmq/websockets for the bridge)
RUN pip3 install --no-cache-dir --break-system-packages numpy scipy pyzmq websockets
ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
WORKDIR /work
COPY docker/entrypoint.sh /usr/local/bin/swash-entrypoint
RUN chmod +x /usr/local/bin/swash-entrypoint
ENTRYPOINT ["/usr/local/bin/swash-entrypoint"]
CMD ["bash"]

# ---------------------------------------------------------------- GPU layer
FROM nvcr.io/nvidia/nvhpc:${NVHPC_TAG} AS gpu
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran gcc g++ cmake ninja-build make perl git \
        python3 python3-pip python3-venv \
        ca-certificates curl less vim-tiny \
    && rm -rf /var/lib/apt/lists/*
RUN pip3 install --no-cache-dir --break-system-packages numpy scipy pyzmq websockets
# nvfortran and the HPC SDK's own MPI are on PATH in the base image
ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
WORKDIR /work
COPY docker/entrypoint.sh /usr/local/bin/swash-entrypoint
RUN chmod +x /usr/local/bin/swash-entrypoint
ENTRYPOINT ["/usr/local/bin/swash-entrypoint"]
CMD ["bash"]
