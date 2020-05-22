//
// Created by diegorenner on 5/21/20.
//

#include "parametrized_mesh.hpp"

#ifndef HELMHOLTZTRANSMISSIONPROBLEM_GEN_SOL_OP_HPP
#define HELMHOLTZTRANSMISSIONPROBLEM_GEN_SOL_OP_HPP

    Eigen::MatrixXcd gen_sol_op(const ParametrizedMesh &mesh,
                               unsigned order,
                               const std::complex<double> k,
                               const double c_o,
                               const double c_i);

    Eigen::MatrixXcd gen_sol_op_1st_der(const ParametrizedMesh &mesh,
                                unsigned order,
                                const std::complex<double> k,
                                const double c_o,
                                const double c_i);

Eigen::MatrixXcd gen_sol_op_2nd_der(const ParametrizedMesh &mesh,
                                    unsigned order,
                                    const std::complex<double> k,
                                    const double c_o,
                                    const double c_i);
#endif //HELMHOLTZTRANSMISSIONPROBLEM_GEN_SOL_OP_HPP
