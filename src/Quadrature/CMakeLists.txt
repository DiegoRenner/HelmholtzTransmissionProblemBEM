include_directories(${PROJECT_SOURCE_DIR}/include)


FILE(GLOB QuadratureSources ./*.cpp)

add_library(quadrature STATIC ${QuadratureSources} ../../include/integral_gauss.hpp)



