cmake_minimum_required(VERSION 3.5)
project(HelmholtzTransmissionProblem)


###########################################################################
#
# INCLUDE MODULES AND FUNCTIONS FOR CMAKE
#
###########################################################################
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
include(${PROJECT_SOURCE_DIR}/cmake/functions/add_all_subdirectories.cmake)

###########################################################################
#
# COMPILATION SETTINGS
#
###########################################################################

# STORE BINARIES IN ./bin
#(else CMake cannot find files to link for automatic targets)
###########################################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# SET DEFAULT BUILD TYPE
###########################################################################
#if( NOT CMAKE_BUILD_TYPE )
#    message(STATUS "No build type selected. Default to Release")
    set(CMAKE_BUILD_TYPE "Release")
#endif()

# SET COMMAND LINE OUPUT IF REQUESTED
###########################################################################
option(CMDL "Command line output requested." OFF)
if (CMDL)
    add_definitions(-DCMDL)
endif()
unset(CMDL CACHE)


# COMPILATION FLAGS
###########################################################################
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

# IGNORE WARNINGS
###########################################################################

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_definitions(-Wno-deprecated-register -Wno-return-type-c-linkage -Wno-keyword-macro -Wno-missing-braces)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    add_definitions(-Wno-deprecated-declarations -Wno-ignored-attributes -Wno-misleading-indentation -Wno-unknown-pragmas)
endif()

###########################################################################
#
# FIND DEPENDENCIES
#
###########################################################################

# FIND EIGEN & COMPLEX BESSEL LIBRARY
###########################################################################

include(${PROJECT_SOURCE_DIR}/cmake/functions/GetEIGEN3.cmake)
add_definitions(-DEIGEN3_ACTIVATED)
include(${PROJECT_SOURCE_DIR}/cmake/functions/Getcomplex_bessel.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/functions/Getarpackpp.cmake)

# ADD SUBDIRECTORIES
add_subdirectory(doxygen)
add_subdirectory(examples)
add_subdirectory(source)
add_subdirectory(test)
