include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(solvers STATIC solvers.cpp)
target_link_libraries(solvers quadrature sol_op)
add_library(compute_SV_der STATIC singular_values.cpp)
target_link_libraries(compute_SV_der quadrature)
add_library(find_zeros STATIC roots.cpp)
add_library(gen_sol_op STATIC gen_sol_op.cpp)
target_link_libraries(gen_sol_op sol_op sol_op_der sol_op_der2)
