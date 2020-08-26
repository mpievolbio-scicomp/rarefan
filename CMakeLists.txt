include(ExternalProject)
include(CMakeDependentOption)
find_program(MAKE_EXE NAMES make gmake)

cmake_minimum_required (VERSION 3.12)

set(CMAKE_DISABLE_SOURCE_CHANGES OFF)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (repinpop)

# Check and locate GSL
set (GSL_ROOT_DIR ${CONDA_PREFIX})
find_package(GSL REQUIRED)

# External codes
# divsufsort
ExternalProject_Add(divsufsort
    GIT_REPOSITORY https://github.com/y-256/libdivsufsort.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    BUILD_IN_SOURCE false
    INSTALL_COMMAND ${MAKE_EXE} install
)

# Count number of processors available for parallel building.
include(ProcessorCount)
ProcessorCount(NPROC)
MESSAGE("Found ${NPROC} processors")
if(NOT NPROC EQUAL 0)
    set(CMAKE_BUILD_FLAGS -j${NPROC})
endif()

# andi
ExternalProject_Add(andi
    URL https://github.com/EvolBioInf/andi/archive/v0.13.tar.gz
    DOWNLOAD_NAME andi-0.13.tar.gz
    CONFIGURE_COMMAND autoreconf -fi -Im4 && ./configure CPPFLAGS=-I$ENV{CONDA_PREFIX}/include CFLAGS=-I$ENV{CONDA_PREFIX}/include LDFLAGS=-L$ENV{CONDA_PREFIX}/lib --prefix=$ENV{CONDA_PREFIX}
    BUILD_IN_SOURCE true
    INSTALL_COMMAND ${MAKE_EXE} install
)

# clustdist
ExternalProject_Add(clustdist
    GIT_REPOSITORY https://github.com/EvolBioInf/clustDist.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND cd ../clustdist && ${MAKE_EXE}
    INSTALL_COMMAND cd ../clustdist && install build/clustDist $ENV{CONDA_PREFIX}/bin
)

# legacy_blast patch
execute_process(patch ${CONDA_PREFIX}/bin/legacy_blast.pl patch)

# provide libRlapack requested by r-ape.
execute_process(ln -s ${CONDA_PREFIX}/lib/R/modules/lapack.so ${CONDA_PREFIX}/lib/libRlapack.so)
