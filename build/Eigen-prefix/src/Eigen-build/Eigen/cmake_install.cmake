# Install script for directory: /home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen_install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Cholesky"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/CholmodSupport"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Core"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Dense"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Eigen"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Eigenvalues"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Geometry"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Householder"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/IterativeLinearSolvers"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Jacobi"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/LU"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/MetisSupport"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/OrderingMethods"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/PaStiXSupport"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/PardisoSupport"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/QR"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/QtAlignedMalloc"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/SPQRSupport"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/SVD"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/Sparse"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/SparseCholesky"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/SparseCore"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/SparseLU"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/SparseQR"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/StdDeque"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/StdList"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/StdVector"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/SuperLUSupport"
    "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/home/luka/Projects/PMF/HelmholtzTransmissionProblemBEM/build/Eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

