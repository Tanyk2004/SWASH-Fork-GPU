# a CMake script to build SWASH
#
# Copyright (C) 2023  Delft University of Technology

# use of CMake version 3.12 or greater
cmake_minimum_required( VERSION 3.12 )

# by default, the SWASH package will be installed in a default location
if( WIN32 )
  set( CMAKE_INSTALL_PREFIX "C:/Program Files/swash" CACHE STRING "" )
else()
  set( CMAKE_INSTALL_PREFIX "/usr/local/swash" CACHE STRING "" )
endif()

# Fortran project for SWASH
project( SWASH Fortran )

# set the SWASH version
set( VERSION 9.01A )

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
