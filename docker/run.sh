#!/usr/bin/env bash
# Run a command inside the SWASH container with the repository mounted at /work.
#
#   docker/run.sh [cpu|gpu] <command...>
#
# Examples
#   docker/run.sh cpu docker/build-all.sh cpu            # gfortran build + goldens
#   docker/run.sh gpu docker/build-all.sh gpu 120        # nvfortran/OpenACC build for an RTX 50xx (cc 120); 89 for RTX 40xx
#   docker/run.sh gpu tests/run_case.sh build-nvhpc marina_lay2 --tier gpu
#   docker/run.sh cpu                                    # interactive shell
#
# Works with docker (rootless is fine) or podman; no sudo needed when your user
# is in the docker group or rootless docker/podman is set up.  GPU access needs
# the NVIDIA Container Toolkit on the host (the admin installs it once).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/.." && pwd)"
target="${1:-cpu}"; shift || true
image="swash:${target}"
engine="${CONTAINER_ENGINE:-}"
if [[ -z "$engine" ]]; then
  if command -v docker >/dev/null 2>&1; then engine=docker; elif command -v podman >/dev/null 2>&1; then engine=podman; else echo "neither docker nor podman found" >&2; exit 1; fi
fi
gpuflags=()
if [[ "$target" == "gpu" ]]; then
  if [[ "$engine" == "podman" ]]; then gpuflags=(--device nvidia.com/gpu=all --security-opt=label=disable); else gpuflags=(--gpus all); fi
fi
tty=(); [[ -t 0 && -t 1 ]] && tty=(-it)
# keep files owned by the invoking user (docker only; podman rootless maps root to the user already)
userflag=(); [[ "$engine" == "docker" ]] && userflag=(--user "$(id -u):$(id -g)" -e HOME=/tmp)
exec "$engine" run --rm "${tty[@]}" "${gpuflags[@]}" "${userflag[@]}" \
  -v "$root":/work -w /work --shm-size=1g "$image" "${@:-bash}"
