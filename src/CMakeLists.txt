# define switches
if( WIN32 )
  set( SWITCHES -dos -cvis )
elseif( UNIX )
  set( SWITCHES -unix )
endif()

# include user-dependent switches
if( MPI )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
    list( APPEND SWITCHES -mpi )
  else()
    list( APPEND SWITCHES -impi -mpi )
  endif()
endif()
if( TIMG )
  # enable the built-in section timers (SWTSTA/SWTSTO/SWPRTI) so that the PRINT file ends with a timing table
  list( APPEND SWITCHES -timg )
endif()
if( OPENACC )
  # enable OpenACC-only code paths (lines prefixed with !ACC); the !$acc directives themselves need no switch
  list( APPEND SWITCHES -acc )
endif()

# set switches
if( PERL_FOUND )
  execute_process( COMMAND ${PERL_EXECUTABLE} ${CMAKE_SOURCE_DIR}/switch.pl ${SWITCHES} ${SRC}/*.ftn ${SRC}/*.ftn90 )
else()
  message( FATAL_ERROR "-- Cannot build ${CMAKE_PROJECT_NAME} without Perl" )
endif()

# set compiler flags
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  # GNU Fortran
  set( CMAKE_Fortran_FLAGS "-O -w -fno-second-underscore -ffree-line-length-none" CACHE STRING "" FORCE )
  if( OPENACC )
    message( WARNING "-- OPENACC=ON with gfortran: directives are compiled for the host only (no GPU offload); use nvfortran for GPU builds" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenacc" )
  endif()
  if( ${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10 )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch" )
  endif()
  if( WIN32 AND ${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 7 )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdec" )
  endif()
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  # Intel Fortran
  if( UNIX )
    set( CMAKE_Fortran_FLAGS "-O2 -W0 -assume byterecl -traceback -diag-disable 8290 -diag-disable 8291 -diag-disable 8293" CACHE STRING "" FORCE )
  elseif( WIN32 )
    set( CMAKE_Fortran_FLAGS "/O2 /assume:byterecl /traceback /nowarn /nologo /Qdiag-disable:8290 /Qdiag-disable:8291 /Qdiag-disable:8293" CACHE STRING "" FORCE )
  endif()
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" )
  # Portland Group / NVIDIA HPC SDK (nvfortran); the only compiler with OpenACC GPU offload support here
  set( CMAKE_Fortran_FLAGS "-O2 -Mbackslash" CACHE STRING "" FORCE )
  if( OPENACC )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc=gpu -gpu=cc${GPU_CC} -Minfo=accel" )
  endif()
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Fujitsu" )
  # Lahey Fortran
  set( CMAKE_Fortran_FLAGS "-O1 -staticlink -nwo" CACHE STRING "" FORCE )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  # IBM XL Fortran
  set( CMAKE_Fortran_FLAGS "-O3 -qstrict -qarch=auto -qtune=auto -qcache=auto -qalign=4k -w" CACHE STRING "" FORCE )
else()
  message( FATAL_ERROR "-- Current Fortran compiler ${CMAKE_Fortran_COMPILER} not supported" )
endif()

# if desired, include MPI parallelisation
if( MPI )
  find_package( MPI REQUIRED COMPONENTS Fortran )
else()
  # turn off MPI
  unset( MPI_Fortran_FOUND CACHE )
  unset( MPI_Fortran_COMPILER CACHE )
  unset( MPI_Fortran_LIBRARIES CACHE )
endif()

# build library
if( NOT TARGET swash${VERSION} )
  include( ${SRC}/srclist.cmake )
  add_library( swash${VERSION} STATIC ${swash_src} )
  set_target_properties( swash${VERSION} PROPERTIES LINKER_LANGUAGE Fortran )
  if( ENGINE )
    # position independent code is required because the engine shared library links this archive
    set_target_properties( swash${VERSION} PROPERTIES POSITION_INDEPENDENT_CODE ON )
  endif()
endif()

# build the engine shared library (C-callable control API around the solver, see engine/)
if( ENGINE AND NOT TARGET swash_engine )
  set( ENGINE_SRC ${CMAKE_SOURCE_DIR}/engine )
  execute_process( COMMAND ${PERL_EXECUTABLE} ${CMAKE_SOURCE_DIR}/switch.pl ${SWITCHES} ${ENGINE_SRC}/*.ftn90 )
  add_library( swash_engine SHARED ${ENGINE_SRC}/swash_state.f90 ${ENGINE_SRC}/swash_ctl.f90 )
  set_target_properties( swash_engine PROPERTIES LINKER_LANGUAGE Fortran POSITION_INDEPENDENT_CODE ON )
  target_link_libraries( swash_engine swash${VERSION} ${MPI_Fortran_LIBRARIES} )
  install( TARGETS swash_engine LIBRARY DESTINATION lib )
  install( FILES ${ENGINE_SRC}/swash_engine.h DESTINATION include )
endif()

# build executable
if( NOT TARGET ${EXE} )
  set( swashmain ${SRC}/Swash.f90 )
  add_executable( ${EXE} ${swashmain} )
  target_link_libraries( ${EXE} swash${VERSION} ${MPI_Fortran_LIBRARIES} )
endif()

# install the package
install( TARGETS ${EXE} swash${VERSION} RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib )
install( DIRECTORY ${CMAKE_BINARY_DIR}/mod/ DESTINATION mod )
install( DIRECTORY ${CMAKE_SOURCE_DIR}/doc/ DESTINATION doc )
install( DIRECTORY ${CMAKE_SOURCE_DIR}/misc/ DESTINATION misc )
install( DIRECTORY ${CMAKE_SOURCE_DIR}/tools/ DESTINATION tools )
if( WIN32 )
  install( PROGRAMS ${CMAKE_SOURCE_DIR}/bin/swashrun.bat DESTINATION bin )
elseif( UNIX )
  install( PROGRAMS ${CMAKE_SOURCE_DIR}/bin/swashrun DESTINATION bin )
endif()
