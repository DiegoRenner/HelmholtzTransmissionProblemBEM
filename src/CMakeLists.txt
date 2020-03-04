set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra")

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(Parametrizations)

add_subdirectory(Quadrature)

add_library(single_layer STATIC single_layer.cpp Parametrizations/parametrized_mesh.cpp)
add_library(double_layer STATIC double_layer.cpp Parametrizations/parametrized_mesh.cpp)
add_library(hypersingular STATIC hypersingular.cpp Parametrizations/parametrized_mesh.cpp)
add_library(adj_double_layer STATIC adj_double_layer.cpp Parametrizations/parametrized_mesh.cpp)
add_library(mass_matrix STATIC mass_matrix.cpp Parametrizations/parametrized_mesh.cpp)
add_library(solvers STATIC solvers.cpp single_layer double_layer hypersingular adj_double_layer mass_matrix)