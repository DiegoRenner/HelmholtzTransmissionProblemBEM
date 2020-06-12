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

#ifndef SINGULAR_VALUESHPP
#define SINGULAR_VALUESHPP

#include "parametrized_mesh.hpp"

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
                   const double *list,
                   const unsigned count);

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
                           const double *list,
                           const unsigned count);
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
                           const double *list,
                           const unsigned count);

/**
 * Extrapolation based numerical differentation
 * with a posteriori error control.
 * This matrix was taken from the lecture document of
 * Prof. R. Hiptmair's course "Advanced Numerical Methods
 * for Computational Science and Engineering"
 * @param f handle of a function defined in a neighbourhood of x ∈ R
 * @param x point at which approximate derivative is desired
 * @param h0 initial distance from x
 * @param rtol relative target tolerance
 * @param atol absolute tolerance
 * @return value of computed derivative
 */
double der_by_ext( std::function<double(double)> f ,
                   const double x ,
                   const double h0 ,
                   const double rtol ,
                   const double atol );

#endif // SINGULAR_VALUESHPP
