/**
 * \file compute.hpp
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem of the form given in \f$\eqref{eq:dirbvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef COMPUTESVDERHPP
#define COMPUTESVDERHPP

#include "parametrized_mesh.hpp"

    Eigen::VectorXd sv(const Eigen::MatrixXcd &T,
                               const double *list,
                               const unsigned count);

    Eigen::MatrixXd sv_1st_der(const Eigen::MatrixXcd &T,
            const Eigen::MatrixXcd &T_der,
            const double *list,
            const unsigned count);
    Eigen::MatrixXd sv_2nd_der(const Eigen::MatrixXcd &T,
                               const Eigen::MatrixXcd &T_der,
                               const Eigen::MatrixXcd &T_der_2,
                               const double *list,
                               const unsigned count);

    // Extrapolation based numerical differentation
    // with a posteriori error control
    // f: handle of a function defined in a neighbourhood of x ∈ R
    // x: point at which approximate derivative is desired
    // h0: initial distance from x
    // rtol: relative target tolerance, atol: absolute tolerance
    double der_by_ext( std::function<double(double)> f ,
                       const double x ,
                       const double h0 ,
                       const double rtol ,
                       const double atol );

#endif // DIRICHLETHPP
