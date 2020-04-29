cmake_minimum_required(VERSION 3.15)
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
if( NOT CMAKE_BUILD_TYPE )
    message(STATUS "No build type selected. Default to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# COMPILATION FLAGS
###########################################################################
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

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

# FIND EIGEN
###########################################################################

include(${PROJECT_SOURCE_DIR}/cmake/functions/GetEIGEN3.cmake)
add_definitions(-DEIGEN3_ACTIVATED)

find_package(Boost REQUIRED COMPONENTS program_options)
if (${Boost_MAJOR_VERSION} VERSION_LESS 1 OR ${Boost_MINOR_VERSION} VERSION_LESS 58)
    message(STATUS "Downloading odeint, since boost version too old.")
    ExternalProject_Add(
            ODEINT
            URL https://github.com/headmyshoulder/odeint-v2/archive/v2.6_boost1.60.zip
            SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/odeint_source
            INSTALL_COMMAND  ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
    )
    include_directories (${CMAKE_CURRENT_BINARY_DIR}/odeint_source/include)
else()
    add_custom_target(ODEINT) # dependency dummy
endif()
# Ugly fix errors with boost versions lower then 1.66
if(${Boost_MAJOR_VERSION} VERSION_LESS 1 OR ${Boost_MINOR_VERSION} VERSION_LESS 66)
    message("Boost version to old to contain chebyshev interpolation. Downloading manually")
    file(DOWNLOAD "https://raw.githubusercontent.com/boostorg/math/develop/include/boost/math/special_functions/chebyshev.hpp" "boost_missing/boost/math/special_functions/chebyshev.hpp")
    include_directories("${CMAKE_CURRENT_BINARY_DIR}/boost_missing")
endif()
include_directories (${Boost_INCLUDE_DIRS})
message(${Boost_INCLUDE_DIRS})

add_subdirectory(doxygen)
add_subdirectory(examples)
add_subdirectory(src)
add_subdirectory(test)
