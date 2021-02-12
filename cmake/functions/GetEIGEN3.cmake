include(ExternalProject) # To download Eigen if it is missing

find_package(EIGEN3)

if(EIGEN3_VERSION VERSION_LESS 3.3.0)
  message("Eigen version 3.3 required, but ${EIGEN3_VERSION} found. Downloading 3.3.4 instead.")
  set(EIGEN3_FOUND FALSE)
endif()


if(${EIGEN3_FOUND})

  include_directories(${EIGEN3_INCLUDE_DIR})
  add_custom_target(Eigen) # dependency dummy

else()
  #  if not found system wide download.
    SET(DOWNLOADING_EIGEN ON)
    message("-- Downloading Eigen3")
    ExternalProject_Add(
        Eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.3.4/eigen-3.3.4.zip
         SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
        DOWNLOAD_NO_PROGRESS 1
        CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Eigen_install)
        include_directories(${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3)

endif()

