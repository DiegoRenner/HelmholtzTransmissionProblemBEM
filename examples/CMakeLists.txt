include_directories(${PROJECT_SOURCE_DIR}/include)

set(CMAKE_CXX_FLAGS "-std=c++11")

add_executable(helmholtz_test helmholtz_circle_test_1.cpp)
add_library(generate_solution STATIC generate_solution.cpp )
add_executable(matlab_read_test matlab_read_test.cpp)

target_link_libraries(helmholtz_test single_layer mass_matrix generate_solution parametrizations quadrature hypersingular double_layer)
