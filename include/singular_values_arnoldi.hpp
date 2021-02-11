/**
 * \file singular_values.hpp
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem
 *
 * This File is a part of the 2D-Parametric BEM package
 * It incorporates an extrpolation function taken from the lectur document
 * "Advanced Numerical Methods for Computational Science and Engineering"
 * by Prof. R. Hiptmair.
 */

#ifndef SINGULAR_VALUES_ARNOLDIHPP
#define SINGULAR_VALUES_ARNOLDIHPP

#include "parametrized_mesh.hpp"
#include "st_vec_storage.hpp"

namespace arnoldi {
/**
 * Compute singular values defined by \p list and \p count of the
 * complex matrix \p T.
 * @param T matrix of which to compute singular values
 * @param list contains index of singular values to compute,
 * 0 is smallest
 * @param count number of singular values to be computed
 * @return singular values of complex matrix \p T defined
 * by \p list and \p count
 */
    Eigen::VectorXd sv(const Eigen::MatrixXcd &T,
                       const unsigned count,
                       const double acc=1e-16);

/**
 * Compute derivative of singular values defined by \p list and \p count of the
 * complex matrix \p T.
 * @param T matrix of which to compute singular values
 * @param T_der derivative of \p T, 0 is smallest
 * @param list index of singular values to be computed
 * @param count number of singular values to be computed
 * @return singular values and their first derivative
 * 1. column singular values
 * 2. column derivatives
 */
    Eigen::MatrixXd sv_1st_der(const Eigen::MatrixXcd &T,
                               const Eigen::MatrixXcd &T_der,
                               const unsigned count,
                               const double acc=1e-16);

/**
 * Compute second derivative of singular values defined by \p list and \p count of the
 * complex matrix \p T.
 * @param T matrix of which to compute singular values
 * @param T_der derivative of \p T
 * @param T_der2 second derivative of \p T
 * @param list index of singular values to be computed, 0 is smallest
 * @param count number of singular values to be computed
 * @return singular values and their first two derivatives
 * 1. column, singular values,
 * 2. column, derivatives,
 * 3. column, second derivatives
 */
    Eigen::MatrixXd sv_2nd_der(const Eigen::MatrixXcd &T,
                               const Eigen::MatrixXcd &T_der,
                               const Eigen::MatrixXcd &T_der2,
                               const unsigned count,
                               const double acc=1e-16);

}
#endif // SINGULAR_VALUESHPP
