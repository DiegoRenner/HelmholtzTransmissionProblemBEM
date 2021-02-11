/**
 * \file gen_sol_op.hpp
 * \brief This file contains functions that compute the approximation
 * of the operator and it's first two derivatives of the second-kind direct BIEs
 * for the Helmholtz transmission problemi using Galerkin BEM.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library
 */
#include "parametrized_mesh.hpp"
#include "arpack/rcompsol.h"

#ifndef ARPP_EIG_INTERFACEHPP
#define ARPP_EIG_INTERFACEHPP

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
void arpp_to_eig(std::complex<double>* in, Eigen::VectorXcd& out);
void eig_to_arpp(Eigen::VectorXcd& in, std::complex<double>* out);
void arpp_to_eig(ARrcCompStdEig<double>& in, Eigen::VectorXd& out_vals, Eigen::MatrixXcd& out_vectors);

#endif
