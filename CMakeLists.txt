include(ExternalProject)
include(CMakeDependentOption)
find_program(MAKE_EXE NAMES make gmake)

cmake_minimum_required (VERSION 3.12)

set(CMAKE_DISABLE_SOURCE_CHANGES OFF)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (repinpop)

# Check and locate GSL
set (GSL_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
find_package(GSL REQUIRED)


# Count number of processors available for parallel building.
include(ProcessorCount)
ProcessorCount(NPROC)
MESSAGE("Found ${NPROC} processors")
if(NOT NPROC EQUAL 0)
    set(CMAKE_BUILD_FLAGS -j${NPROC})
endif()

# clustdist
ExternalProject_Add(clustdist
    GIT_REPOSITORY https://github.com/EvolBioInf/clustDist.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND cd ../clustdist && ${MAKE_EXE}
    INSTALL_COMMAND cd ../clustdist && install build/clustDist ${CMAKE_INSTALL_PREFIX}/bin
)

# legacy_blast patch
add_custom_target(patch_blast ALL
	COMMAND patch ${CMAKE_INSTALL_PREFIX}/bin/legacy_blast.pl ../patch
)


# provide libRlapack requested by r-ape.
add_custom_target(link_lapack ALL
	COMMAND ln -sf ${CMAKE_INSTALL_PREFIX}/lib/R/modules/lapack.so ${CMAKE_INSTALL_PREFIX}/lib/libRlapack.so
)

# java code
ExternalProject_Add(repin_ecology
  DOWNLOAD_COMMAND cp -ar ${PROJECT_SOURCE_DIR}/REPIN_ecology/REPIN_ecology/ . && pwd
  CONFIGURE_COMMAND ""
  BUILD_COMMAND cd ../REPIN_ecology && gradle build
  INSTALL_COMMAND install ../REPIN_ecology/build/libs/REPIN_ecology.jar ${CMAKE_INSTALL_PREFIX}/lib/REPIN_ecology.jar
  )






