cmake_minimum_required(VERSION 3.5)

project(opensc)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
 set(dialect "-Wuninitialized -fbacktrace -ffree-form -ffree-line-length-none -std=f2008 -fimplicit-none")
   #     set(dialect "-Wconversion -ffree-form -ffree-line-length-none   -std=f2008 -fimplicit-none")
 #set(dialect "-ffree-form -ffree-line-length-none   -std=f2008 -fimplicit-none")   
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

enable_language(C)
#if(CMAKE_C_COMPILER_ID MATCHES "GNU")
# set(dialect "-Wall ")
#endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wuninitialized -pedantic -Wextra")


#set(CMAKE_BUILD_TYPE Debug)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_Fortran_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_Fortran_FLAGS}")
#endif()


#--------------------------
option(GetFFTW "GetFFTW" ON)
if(GetFFTW)

# FFTW2 
include(ExternalProject)
ExternalProject_Add(project_fftw
  URL http://www.fftw.org/fftw-2.1.5.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
  CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/fftw/src/project_fftw/configure --enable-mpi --prefix=<INSTALL_DIR>
  BUILD_COMMAND make -j 8
  INSTALL_COMMAND make install
  PREFIX=${CMAKE_CURRENT_BINARY_DIR}/fftw
)

ExternalProject_Get_Property(project_fftw install_dir)
message(STATUS "fftw install_dir: " ${install_dir})
add_library(fftw STATIC IMPORTED)
set_property(TARGET fftw PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfftw.a)

add_library(fftw_mpi STATIC IMPORTED)
set_property(TARGET fftw_mpi PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfftw_mpi.a)

add_dependencies(fftw project_fftw)
include_directories(${install_dir}/include)


endif(GetFFTW)
#--------------------------


##find_library (FFTW_LIBRARIES NAMES libdfftw.a)
##message(STATUS "FFTW_LIBRARIES: " ${FFTW_LIBRARIES})

##find_path (FFTW_INCLUDES NAMES fftw)
##message(STATUS "FFTW_INCLUDES: " ${FFTW_INCLUDES})



find_package(MPI REQUIRED)


#--------------------------
#
# Preprocessor
#
option(MPIPARALLEL "use MPI" ON)
IF(MPIPARALLEL)
    ADD_DEFINITIONS(-DMPIPARALLEL)
ENDIF(MPIPARALLEL)

option(FFT_FFTW "use fftw" ON)
IF(FFT_FFTW)
    ADD_DEFINITIONS(-DFFT_FFTW)
ENDIF(FFT_FFTW)

#set (CMAKE_Fortran_FLAGS "${CMAKE_C_FLAGS} ${MPI_Fortran_FLAGS}")
message(STATUS "MPI_Fortran_INCLUDE_PATH: " ${MPI_Fortran_INCLUDE_PATH})

include_directories(SYSTEM ${MPI_Fortran_INCLUDE_PATH})


#--------------------------
#
# Compile.
#
file(GLOB_RECURSE sources  code/*.f90 code/fft/*.h  code/fft/*.c )
add_executable(test_opensc ${sources})
target_link_libraries(test_opensc fftw fftw_mpi ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} )

# Link FFTW
## target_link_libraries(test_opensc ${FFTW_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} )
