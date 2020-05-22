set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra")

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(Parametrizations)

add_subdirectory(Quadrature)

add_library(single_layer STATIC single_layer.cpp Parametrizations/parametrized_mesh.cpp)
add_library(single_layer_der STATIC single_layer_der.cpp Parametrizations/parametrized_mesh.cpp)
add_library(single_layer_der2 STATIC single_layer_der2.cpp Parametrizations/parametrized_mesh.cpp)

add_library(double_layer STATIC double_layer.cpp Parametrizations/parametrized_mesh.cpp)
add_library(double_layer_der STATIC double_layer_der.cpp Parametrizations/parametrized_mesh.cpp)
add_library(double_layer_der2 STATIC double_layer_der2.cpp Parametrizations/parametrized_mesh.cpp)

add_library(hypersingular STATIC hypersingular.cpp Parametrizations/parametrized_mesh.cpp)
add_library(hypersingular_der STATIC hypersingular_der.cpp Parametrizations/parametrized_mesh.cpp)
add_library(hypersingular_der2 STATIC hypersingular_der2.cpp Parametrizations/parametrized_mesh.cpp)

add_library(adj_double_layer STATIC adj_double_layer.cpp Parametrizations/parametrized_mesh.cpp)

add_library(mass_matrix STATIC mass_matrix.cpp Parametrizations/parametrized_mesh.cpp)

add_library(solvers STATIC Utilities/solvers.cpp Quadrature/logweight_quadrature.cpp Quadrature/genLaguerreRule.cpp
        single_layer
        double_layer adj_double_layer hypersingular mass_matrix)
add_library(compute_SV_der STATIC Utilities/singular_values.cpp Quadrature/logweight_quadrature.cpp Quadrature/genLaguerreRule.cpp
        single_layer single_layer_der single_layer_der2
        double_layer double_layer_der double_layer_der2
        hypersingular hypersingular_der hypersingular_der2
        adj_double_layer mass_matrix)
add_library(find_zeros STATIC Utilities/roots.cpp)
