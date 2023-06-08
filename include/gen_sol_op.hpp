/**
 * \file gen_sol_op.hpp
 * \brief This file contains functions that compute the approximation
 * of the operator and it's first two derivatives of the second-kind direct BIEs
 * for the Helmholtz transmission problemi using Galerkin BEM.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library
 */

#ifndef GEN_SOL_OPHPP
#define GEN_SOL_OPHPP

#include "parametrized_mesh.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "gauleg.hpp"

class SolutionsOperator
{
private:
    Eigen::VectorXcd L, U;
    unsigned int numpanels, order;
    DiscontinuousSpace<0> discont_space;
    ContinuousSpace<1> cont_space;
    const ParametrizedMesh &mesh;
    Eigen::MatrixXcd M;
    QuadRule GaussQR;
    QuadRule CGaussQR;
    /**
     * LDLT decomposition of tridiagonal matrix A.
     * Lower subdiagonal is stored to l, diagonal do d.
     */
    void tridiagonal_ldl(const Eigen::MatrixXcd &A, Eigen::VectorXcd &l, Eigen::VectorXcd &d) const;
    /**
     * Solve A*X = B where A is lower triangular Hessenberg
     * and L, U is its LU decomposition.
     */
    Eigen::MatrixXcd tridiagonal_lu(const Eigen::MatrixXcd &B) const;

public:
    /**
     * Initialize solutions operator class.
     * @param mesh boundary on which to compute transmission problem
     * @param order order of quadrature rule with which to compute
     * Galerking matrix entries
     */
    SolutionsOperator(const ParametrizedMesh &mesh_in, unsigned order_in);
    /**
     * Compute approximation of solutions operator for second-kind direct BIEs of
     * Helmholtz transmission problem using Galerkin BEM.
     * @param k wavenumber
     * @param c_o refraction index of outer domain
     * @param c_i refraction indef of inner domain
     * @return solutions operator approximation
     */
    Eigen::MatrixXcd gen_sol_op(const std::complex<double> &k, double c_o, double c_i) const;
    /**
     * Compute approximation of derivative of
     * solutions operator for second-kind direct BIEs of
     * Helmholtz transmission problem using Galerkin BEM.
     * @param k wavenumber
     * @param c_o refraction index of outer domain
     * @param c_i refraction indef of inner domain
     * @return derivative solutions operator approximation
     */
    Eigen::MatrixXcd gen_sol_op_1st_der(const std::complex<double> &k, double c_o, double c_i) const;
    /**
     * Compute approximation of second derivative of
     * solutions operator for second-kind direct BIEs of
     * Helmholtz transmission problem using Galerkin BEM.
     * @param k wavenumber
     * @param c_o refraction index of outer domain
     * @param c_i refraction indef of inner domain
     * @return second derivative of solutions operator approximation
     */
    Eigen::MatrixXcd gen_sol_op_2nd_der(const std::complex<double> &k, double c_o, double c_i) const;
    /**
     * The following routines construct matrices K, V and W separately.
     */
    Eigen::MatrixXcd K_cont(const std::complex<double> &k, double c) const;
    Eigen::MatrixXcd K_cont_discont(const std::complex<double> &k, double c) const;
    Eigen::MatrixXcd K_discont_cont(const std::complex<double> &k, double c) const;
    Eigen::MatrixXcd V_cont(const std::complex<double> &k, double c) const;
    Eigen::MatrixXcd V_discont(const std::complex<double> &k, double c) const;
    Eigen::MatrixXcd W_cont(const std::complex<double> &k, double c) const;
    /**
     * Return mass matrix
     */
    const Eigen::MatrixXcd &mass_matrix_cont() const { return M; };
};

#endif //GEN_SOL_OPHPP
