/**
 * \file dirichlet.hpp
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem of the form given in \f$\eqref{eq:dirbvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef TRANSMISSIONHPP
#define TRANSMISSIONHPP

#include "parametrized_mesh.hpp"

namespace parametricbem2d {
    namespace bvp {
/**
 * This namespace contains the solver using the direct first kind method which
 * has the variational formulation as given in \f$\eqref{eq:aVdir}\f$. The
 * Solver uses the lowest order BEM spaces for computation.
 */
        namespace direct_first_kind {
            Eigen::VectorXcd solve_dirichlet(const ParametrizedMesh &mesh,
                                             std::function<std::complex<double>(double, double)> u_dir,
                                             std::function<std::complex<double>(double, double)> u_neu,
                                             unsigned order,
                                             const double k);

            Eigen::VectorXcd solve_neumann(const ParametrizedMesh &mesh,
                                           std::function<std::complex<double>(double, double)> u_dir,
                                           std::function<std::complex<double>(double, double)> u_neu,
                                           unsigned order,
                                           const double k);
        } // namespace direct_first_kind
    } // namespace bvp
    namespace tsp {
        namespace direct_second_kind {
            Eigen::VectorXcd solve_debug_3(const ParametrizedMesh &mesh,
                                           std::function<std::complex<double>(double, double)> u_inc_dir,
                                           std::function<std::complex<double>(double, double)> u_inc_neu,
                                           std::function<std::complex<double>(double, double)> sol_dir,
                                           std::function<std::complex<double>(double, double)> sol_neu,
                                           unsigned order,
                                           const double k_o,
                                           const double k_i);
            Eigen::VectorXcd solve_debug_2(const ParametrizedMesh &mesh,
                                           std::function<std::complex<double>(double, double)> u_inc_dir,
                                           std::function<std::complex<double>(double, double)> u_inc_neu,
                                           std::function<std::complex<double>(double, double)> sol_dir,
                                           std::function<std::complex<double>(double, double)> sol_neu,
                                           unsigned order,
                                           const double k);

            Eigen::VectorXcd solve_debug_1(const ParametrizedMesh &mesh,
                                           std::function<std::complex<double>(double, double)> u_inc_dir,
                                           std::function<std::complex<double>(double, double)> u_inc_neu,
                                           std::function<std::complex<double>(double, double)> sol_dir,
                                           std::function<std::complex<double>(double, double)> sol_neu,
                                           unsigned order,
                                           const double k);

            Eigen::VectorXcd solve(const ParametrizedMesh &mesh,
                                   std::function<std::complex<double>(double, double)> u_inc_dir,
                                   std::function<std::complex<double>(double, double)> u_inc_neu,
                                   std::function<std::complex<double>(double, double)> sol_dir,
                                   std::function<std::complex<double>(double, double)> sol_neu,
                                   unsigned order,
                                   const double k_o,
                                   const double k_i);
            Eigen::MatrixXcd compute_operator(const ParametrizedMesh &mesh,
                                   unsigned order,
                                   std::complex<double> k_o,
                                   std::complex<double> k_i);
        } // namespace direct_second_kind
/**
 * This function is used to solve the Dirichlet boundary value problem given
 * in \f$\eqref{eq:dirbvp}\f$ using the variational formulation given in
 * \f$\eqref{eq:l2dv}\f$ after lifting of the variation formulation given in
 * \f$\eqref{eq:bie2nbpvv}\f$ to the \f$L^{2}(\Gamma)\f$ space. The function
 * outputs a vector of estimated Neumann trace of the solution u.
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param g Dirichlet boundary condition in 2D using a function of the form
 *          double(double,double)
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd type representing the Neumann trace of the
 * solution u
 */
    } // namespace tsp
/**
 * This namespace contains the solver using the indirect first kind method which
 * has the variational formulation as given in \f$\eqref{eq:iddirVv}\f$. The
 * Solver uses the lowest order BEM spaces for computation.
 */
} // namespace parametricbem2d
#endif // DIRICHLETHPP
