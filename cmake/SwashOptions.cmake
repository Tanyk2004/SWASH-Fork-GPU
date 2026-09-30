# Additional build options for the GPU / engine fork of SWASH
#
#   OPENACC  - compile the OpenACC directives for GPU offload (needs nvfortran from the NVIDIA HPC SDK)
#   GPU_CC   - GPU compute capability passed to -gpu=cc<GPU_CC> (e.g. 80 A100, 89 RTX 4090, 90 H100, 120 RTX 5090)
#   TIMG     - switch on SWASH's built-in section timers; the PRINT file then ends with a timing table
#   ENGINE   - also build libswash_engine (shared library with a C control API, see engine/)
#
# Typical configurations:
#   cmake .. -G Ninja                                           CPU reference build (gfortran)
#   cmake .. -G Ninja -DTIMG=ON                                 CPU build with timers (profiling baseline)
#   cmake .. -G Ninja -DCMAKE_Fortran_COMPILER=nvfortran \
#            -DOPENACC=ON -DGPU_CC=89 -DTIMG=ON                 GPU build
#   cmake .. -G Ninja -DENGINE=ON                               CPU build + engine shared library

option( OPENACC "enable OpenACC GPU offload (nvfortran)" OFF )
option( TIMG    "enable built-in section timers"          OFF )
option( ENGINE  "build the libswash_engine shared library" OFF )
set( GPU_CC "native" CACHE STRING "GPU compute capability for -gpu=cc<value> (e.g. 80, 86, 89, 90, 120, or native)" )

if( OPENACC AND MPI )
  message( STATUS "-- OPENACC with MPI: halo exchanges are host-staged (see docs/BUILD_GPU.md)" )
endif()
