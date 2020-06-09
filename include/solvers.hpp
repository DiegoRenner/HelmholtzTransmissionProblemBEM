/**
 * \file
 * \brief This file defines lowest order solvers for direct first kind BIEs
 * for the Dirichlet and the Neumann problem as well as a direct second kind
 * BIE for the Helmholtz transmission problem.
 *
 * This File is a part of the HelmholtzTransmissionBEM
 */

#ifndef SOLVERSHPP
#define SOLVERSHPP

#include "parametrized_mesh.hpp"
/**
 * This namespace contains solvers for boundary value problems
 */
namespace bvp {

/**
 * This namespace contains solvers using direct first kind BIE
 * for the Dirichlet and the Neumann problem.
 * The solvers use the lowest order BEM spaces for computation.
 */
    namespace direct_first_kind {

        /**
         * This function solves the Dirichlet problem given a mesh \p mesh, the Dirichlet data of u
         * \p u_dir, the order of the quadrature ruile used to compute the Galerkin matrix entries
         * \p order and the wavenumber \p k.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param u_dir Dirichlet data
         * @param order order of qudrature rule for matrix entries
         * @param k wavenumber
         * @return Neumann data of u
         */
        Eigen::VectorXcd solve_dirichlet(const ParametrizedMesh &mesh,
                                         const std::function<std::complex<double>(double, double)> u_dir,
                                         const unsigned order,
                                         const double k);

        /**
         * This function solves the Neumann problem given a mesh \p mesh, the Neumann data of u
         * \p u_dir, the order of the quadrature ruile used to compute the Galerkin matrix entries
         * \p order and the wavenumber \p k.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param u_dir Dirichlet data
         * @param order order of qudrature rule for matrix entries
         * @param k wavenumber
         * @return Dirichlet data of u
         */
        Eigen::VectorXcd solve_neumann(const ParametrizedMesh &mesh,
                                       const std::function<std::complex<double>(double, double)> u_neu,
                                       const unsigned order,
                                       const double k);
    } // namespace direct_first_kind
} // namespace bvp

/**
 * This namespace contains the solver for our transmission problem
 */
namespace tp {

    /**
     * This namespace contains the solver using direct second kind BIEs
     * for the Helmholtz transmission problem.
     * The solver uses the lowest order BEM spaces for computation.
     */
    namespace direct_second_kind {
        /**
         * This function returns the solution to the Helmholtz transmission problem
         * on boundary given by \p mesh for an incoming wave defined by \u_inc_dir and
         * \p u_inc_neu. The wavenumber is set by \p k and th refraction indeces by
         * \p c_o and \p c_i. The Galerkin matrix entries are compute with a quadrature rule
         * defined by the parameter \p order.
         * @param mesh mesh of the boundary on which to compute BIOs
         * @param u_inc_dir Dirichlet data of incoming wave
         * @param u_inc_neu Neumann data of incoming wave
         * @param order order of qudrature rule for matrix entries
         * @param k wavenumber
         * @param c_o refraction index outer domain
         * @param c_i refraction index on inner domain
         * @return Dirichlet and Neumann data of resulting wave
         */
        Eigen::VectorXcd solve(const ParametrizedMesh &mesh,
                               const std::function<std::complex<double>(double, double)> u_inc_dir,
                               const std::function<std::complex<double>(double, double)> u_inc_neu,
                               const unsigned order,
                               const double k,
                               const double c_o,
                               const double c_i);
    } // namespace direct_second_kind
} // namespace tp
#endif // DIRICHLETHPP
