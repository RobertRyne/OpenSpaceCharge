cmake_minimum_required(VERSION 3.5)

project(hello)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
##    set(dialect "-Wall -ffree-form -ffree-line-length-none -std=f2008 -fimplicit-none")
        set(dialect "-Wconversion -ffree-form -ffree-line-length-none   -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compile.
#
file(GLOB_RECURSE sources  code/*.f90 code/*.h)
add_executable(test_opensc ${sources})
add_library(opensc code/open_spacecharge_mod.f90)
