/**
 * \file gen_sol_op.hpp
 * \brief This file contains functions that compute the approximation
 * of the operator and it's first two derivatives of the second-kind direct BIEs
 * for the Helmholtz transmission problemi using Galerkin BEM.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library
 */
#include "parametrized_mesh.hpp"

#ifndef GEN_SOL_OPHPP
#define GEN_SOL_OPHPP

/**
 * Compute approximation of solutions operator for second-kind direct BIEs of
 * Helmholtz transmission problem using Galerkin BEM.
 * @param mesh boundary on which to compute transmission problem
 * @param order order of quadrature rule with which to compute
 * Galerking matrix entries
 * @param k wavenumber
 * @param c_o refraction index of outer domain
 * @param c_i refraction indef of inner domain
 * @return solutions operator approximation
 */
Eigen::MatrixXcd gen_sol_op(const ParametrizedMesh &mesh,
                            unsigned order,
                            const std::complex<double> k,
                            const double c_o,
                            const double c_i);

/**
 * Compute approximation of derivative of
 * solutions operator for second-kind direct BIEs of
 * Helmholtz transmission problem using Galerkin BEM.
 * @param mesh boundary on which to compute transmission problem
 * @param order order of quadrature rule with which to compute
 * Galerking matrix entries
 * @param k wavenumber
 * @param c_o refraction index of outer domain
 * @param c_i refraction indef of inner domain
 * @return derivative solutions operator approximation
 */
Eigen::MatrixXcd gen_sol_op_1st_der(const ParametrizedMesh &mesh,
                                    unsigned order,
                                    const std::complex<double> k,
                                    const double c_o,
                                    const double c_i);

/**
 * Compute approximation of second derivative of
 * solutions operator for second-kind direct BIEs of
 * Helmholtz transmission problem using Galerkin BEM.
 * @param mesh boundary on which to compute transmission problem
 * @param order order of quadrature rule with which to compute
 * Galerking matrix entries
 * @param k wavenumber
 * @param c_o refraction index of outer domain
 * @param c_i refraction indef of inner domain
 * @return second derivative of solutions operator approximation
 */
Eigen::MatrixXcd gen_sol_op_2nd_der(const ParametrizedMesh &mesh,
                                    unsigned order,
                                    const std::complex<double> k,
                                    const double c_o,
                                    const double c_i);
#endif //GEN_SOL_OPHPP
