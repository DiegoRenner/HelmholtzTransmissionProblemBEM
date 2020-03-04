include_directories(${PROJECT_SOURCE_DIR}/include)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra")

add_executable(helmholtz_test helmholtz_circle_test.cpp)
add_executable(matlab_read_test matlab_read_test.cpp)
add_executable(dirichlet_example dirichlet_example.cpp)
add_executable(neumann_example neumann_example.cpp)
add_executable(transmission_example transmission_example.cpp)

add_library(generate_solution STATIC generate_solution.cpp )

target_link_libraries(helmholtz_test single_layer mass_matrix generate_solution parametrizations quadrature hypersingular double_layer)
target_link_libraries(dirichlet_example parametrizations generate_solution solvers)
target_link_libraries(neumann_example parametrizations generate_solution solvers)
target_link_libraries(transmission_example parametrizations generate_solution solvers)

