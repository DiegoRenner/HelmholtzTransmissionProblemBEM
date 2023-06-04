include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(solvers STATIC solvers.cpp)
target_link_libraries(solvers quadrature sol_op)

add_library(compute_SV_der STATIC singular_values.cpp)
target_link_libraries(compute_SV_der quadrature)

add_library(compute_SV_der_arnoldi STATIC singular_values_arnoldi.cpp)
target_link_libraries(compute_SV_der_arnoldi arpp_eig_interface)

add_library(gen_sol_op STATIC gen_sol_op.cpp)
target_link_libraries(gen_sol_op sol_op sol_op_der sol_op_der2)


add_library(find_roots STATIC find_roots.cpp)

add_library(arpp_eig_interface STATIC arpp_eig_interface.cpp)

add_library(generate_solution STATIC gen_sol.cpp)

add_library(complex_bessel STATIC cbessel.cpp)

add_library(randomized_svd STATIC randsvd.cpp)
