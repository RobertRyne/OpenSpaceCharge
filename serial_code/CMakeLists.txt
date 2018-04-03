cmake_minimum_required(VERSION 3.5)

project(opensc)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
##    set(dialect "-Wall -ffree-form -ffree-line-length-none -std=f2008 -fimplicit-none")
    set(dialect "-Wconversion -ffree-form -ffree-line-length-none  -cpp -std=f2008 -fimplicit-none") 
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


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

#--------------------------------------------
#--------------------------------------------
option(use_fftw "USE_FFTW" OFF)
if(use_fftw)
    add_definitions(-DUSE_FFTW)
    message(STATUS "Using FFTW")

#--------------------------------------------
option(GetFFTW "GetFFTW" ON)
if(GetFFTW)
# FFTW3 
include(ExternalProject)
ExternalProject_Add(project_fftw
  URL http://www.fftw.org/fftw-3.3.7.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
  CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/project_fftw/configure --enable-openmp --prefix=<INSTALL_DIR>
  BUILD_COMMAND make -j 8
  INSTALL_COMMAND make install
  PREFIX=${CMAKE_CURRENT_BINARY_DIR}/fftw
)

ExternalProject_Get_Property(project_fftw install_dir)
message(STATUS "fftw install_dir: " ${install_dir})

# Libraries
add_library(fftw STATIC IMPORTED)
set_property(TARGET fftw PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfftw3.a)

add_library(fftw_omp STATIC IMPORTED)
set_property(TARGET fftw_omp PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfftw3_omp.a)


add_dependencies(fftw project_fftw)
include_directories(${install_dir}/include)

endif(GetFFTW)
#--------------------------------------------
endif()
#--------------------------------------------
#--------------------------------------------


#
# Compile.
#

#file(GLOB_RECURSE sources  code/*.f90 code/*.h)

set(sources
    code/fast_fourier_am.f90
    code/fft_interface_mod.f90
    code/open_spacecharge_mod.f90
    code/open_spacecharge_mod1.f90
    code/test_mod.f90
    code/test_opensc1.f90
)

add_executable(test_opensc ${sources})

# Link FFTW
#
if(GetFFTW)
    add_dependencies(test_opensc fftw)
    target_link_libraries(test_opensc fftw fftw_omp)
endif()

# Library
add_library(openspacecharge STATIC 
    code/fast_fourier_am.f90
    code/fft_interface_mod.f90
    code/open_spacecharge_mod.f90
    code/open_spacecharge_mod1.f90)
    
        

