#!/usr/bin/env bash
# Entry point: prints the toolchain that is available, then runs the command.
if command -v nvfortran >/dev/null 2>&1; then
  echo "swash container: nvfortran $(nvfortran --version 2>/dev/null | sed -n 2p)"
  if command -v nvidia-smi >/dev/null 2>&1; then nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null | sed 's/^/swash container: GPU /'; else echo "swash container: no GPU visible (run with --gpus all)"; fi
fi
echo "swash container: $(gfortran --version | head -1)"
exec "$@"
