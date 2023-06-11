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

#include "galerkin_matrix_builder.hpp"

class SolutionsOperator
{
    const ParametrizedMesh &mesh;
    const AbstractBEMSpace &test_space;
    const AbstractBEMSpace &trial_space;
    size_t dim_test, dim_trial;
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu;
    Eigen::MatrixXd M; // mass matrix
    QuadRule GaussQR;
    QuadRule CGaussQR;
    Eigen::MatrixXcd project(const Eigen::MatrixXcd &T) const;
    void gen_sol_op_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                       Eigen::MatrixXcd &T) const;
    void gen_sol_op_1st_der_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                               Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der) const;
    void gen_sol_op_2nd_der_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                               Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2) const;

public:
    /**
     * Initialize solutions operator class.
     */
    SolutionsOperator(const ParametrizedMesh &mesh_in,
                      unsigned int order,
                      const AbstractBEMSpace &test_space_in,
                      const AbstractBEMSpace &trial_space_in);
    /**
     * Compute approximation of solutions operator for second-kind direct BIEs of
     * Helmholtz transmission problem using Galerkin BEM.
     * @param k wavenumber
     * @param c_o refraction index of outer domain
     * @param c_i refraction indef of inner domain
     */
    void gen_sol_op(const std::complex<double> &k, double c_o, double c_i,
                    Eigen::MatrixXcd &T) const;
    void gen_sol_op(GalerkinMatrixBuilder &builder, const std::complex<double> &k, double c_o, double c_i,
                    Eigen::MatrixXcd &T) const;
    /**
     * Compute approximation of solutions operator and its 1st derivative
     * for second-kind direct BIEs of Helmholtz transmission problem using Galerkin BEM.
     * @param k wavenumber
     * @param c_o refraction index of outer domain
     * @param c_i refraction indef of inner domain
     */
    void gen_sol_op_1st_der(const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der) const;
    void gen_sol_op_1st_der(GalerkinMatrixBuilder &builder, const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der) const;
    /**
     * Compute approximation of solutions operator and its 1st and 2nd derivative
     * for second-kind direct BIEs of Helmholtz transmission problem using Galerkin BEM.
     * @param k wavenumber
     * @param c_o refraction index of outer domain
     * @param c_i refraction indef of inner domain
     */
    void gen_sol_op_2nd_der(const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2) const;
    void gen_sol_op_2nd_der(GalerkinMatrixBuilder &builder, const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2) const;
    /**
     * Return mass matrix.
     */
    const Eigen::MatrixXd &mass_matrix() const { return M; };
    const QuadRule &get_GaussQR() const { return GaussQR; }
    const QuadRule &get_CGaussQR() const { return CGaussQR; }
};

#endif //GEN_SOL_OPHPP
