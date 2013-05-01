cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  set(CMAKE_Fortran_FLAGS "-g -fno-range-check -mcmodel=large -m64 -funroll-all-loops -fno-f2c -O3")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  set(CMAKE_Fortran_FLAGS "-g -r8 -shared-intel -mcmodel=large -m64")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

file(GLOB F_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.f)

add_library(fortPMF ${F_SOURCES})

install(TARGETS fortPMF DESTINATION lib)
