# a CMake script to build SWASH
#
# Copyright (C) 2026  Delft University of Technology

# use of CMake version 3.20 or greater
cmake_minimum_required( VERSION 3.20 )

# by default, the SWASH package will be installed in user's local or home folder
if( WIN32 )
  cmake_path( SET LOCALSWASH $ENV{LOCALAPPDATA}/Programs/wavemodels/swash )
  set( CMAKE_INSTALL_PREFIX "${LOCALSWASH}" CACHE STRING "" )
else()
  set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/wavemodels/swash" CACHE STRING "" )
endif()

# Fortran project for SWASH
project( SWASH Fortran )

# set the SWASH version
set( VERSION 12.01 )

# enable common build
set( CMAKE_BUILD_TYPE None CACHE STRING "" FORCE )

# set the executable name
if( WIN32 )
   set( EXE swash )
else()
   set( EXE swash.exe )
endif()

# set the Fortran77 extension
if( WIN32 )
   set( EXTF for )
else()
   set( EXTF f   )
endif()

# extra build options for the user
option( MPI "" OFF )
include( ${CMAKE_SOURCE_DIR}/cmake/SwashOptions.cmake )

# define bin/lib//mod/src directories
set( BIN ${CMAKE_BINARY_DIR}/bin )
set( LIB ${CMAKE_BINARY_DIR}/lib )
set( MOD ${CMAKE_BINARY_DIR}/mod )
set( SRC ${CMAKE_SOURCE_DIR}/src )

# have the bin and lib folders placed in the build folder
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB} )

# have the modules placed in the mod folder
set( CMAKE_Fortran_MODULE_DIRECTORY ${MOD} )

# Perl is required for setting switches in the source code
find_package( Perl REQUIRED )

# add the source to the library and place it in the lib folder
add_subdirectory( ${SRC} ${LIB} )

# create the SWASH executable and place it in the bin folder
add_subdirectory( ${SRC} ${BIN} )
