/**
 * \file gen_sol_op.hpp
 * \brief This file contains functions that compute the approximation
 * of the operator and it's first two derivatives of the second-kind direct BIEs
 * for the Helmholtz transmission problemi using Galerkin BEM.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library
 * (c) 2023 Luka Marohnić
 */

#ifndef GEN_SOL_OPHPP
#define GEN_SOL_OPHPP

#include "galerkin_matrix_builder.hpp"

class SolutionsOperator
{
    const BuilderData &builder_data;
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu;
    Eigen::MatrixXd M; // mass matrix
    Eigen::MatrixXcd project(const Eigen::MatrixXcd &T) const;
    size_t dim_test, dim_trial;
    // timing info
    bool profiling;
    unsigned count[3] = {0, 0, 0};
    unsigned total_assembly_time[3] = {0, 0, 0};
    unsigned total_hankel_computation_time[3] = {0, 0, 0};
    unsigned total_interaction_matrix_assembly_time[3] = {0, 0, 0};
    unsigned total_projection_time[3] = {0, 0, 0};
    // solutions operator matrix assembly routines
    void gen_sol_op_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                       Eigen::MatrixXcd &T);
    void gen_sol_op_1st_der_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                               Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der);
    void gen_sol_op_2nd_der_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                               Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2);

public:
    /**
     * Initialize solutions operator class.
     * @param builder_data_in bulder data object
     * @param profiling_in whether to do time profiling
     */
    SolutionsOperator(const BuilderData &builder_data_in, bool profiling_in);
    // destructor
    ~SolutionsOperator();
    /**
     * Compute approximation of solutions operator for second-kind direct BIEs of
     * Helmholtz transmission problem using Galerkin BEM.
     * This routine will create a temporary Galerkin matrix builder.
     * @param k wavenumber (complex)
     * @param c_o refraction index of outer domain (should be 1)
     * @param c_i refraction indef of inner domain (must not be smaller than c_o)
     * @param T complex matrix to which the solutions operator matrix will be stored
     */
    void gen_sol_op(const std::complex<double> &k, double c_o, double c_i,
                    Eigen::MatrixXcd &T);
    /**
     * Compute approximation of solutions operator for second-kind direct BIEs of
     * Helmholtz transmission problem using Galerkin BEM.
     * This routine will use the given Galerkin matrix builder.
     * @param builder Galerkin matrix builder
     * @param k wavenumber (complex)
     * @param c_o refraction index of outer domain (should be 1)
     * @param c_i refraction indef of inner domain (must not be smaller than c_o)
     * @param T complex matrix to which the solutions operator matrix will be stored
     */
    void gen_sol_op(GalerkinMatrixBuilder &builder, const std::complex<double> &k, double c_o, double c_i,
                    Eigen::MatrixXcd &T);
    /**
     * Compute approximation of solutions operator and its 1st derivative
     * for second-kind direct BIEs of Helmholtz transmission problem using Galerkin BEM.
     * This routine will create a temporary Galerkin matrix builder.
     * @param k wavenumber (complex)
     * @param c_o refraction index of outer domain (should be 1)
     * @param c_i refraction indef of inner domain (must not be smaller than c_o)
     * @param T complex matrix to which the solutions operator matrix will be stored
     * @param T_der complex matrix to which the 1st derivative of the solutions operator matrix will be stored
     */
    void gen_sol_op_1st_der(const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der);
    /**
     * Compute approximation of solutions operator and its 1st derivative
     * for second-kind direct BIEs of Helmholtz transmission problem using Galerkin BEM.
     * This routine will use the given Galerkin matrix builder.
     * @param builder Galerkin matrix builder
     * @param k wavenumber (complex)
     * @param c_o refraction index of outer domain (should be 1)
     * @param c_i refraction indef of inner domain (must not be smaller than c_o)
     * @param T complex matrix to which the solutions operator matrix will be stored
     * @param T_der complex matrix to which the 1st derivative of the solutions operator matrix will be stored
     */
    void gen_sol_op_1st_der(GalerkinMatrixBuilder &builder, const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der);
    /**
     * Compute approximation of solutions operator and its 1st and 2nd derivatives
     * for second-kind direct BIEs of Helmholtz transmission problem using Galerkin BEM.
     * This routine will create a temporary Galerkin matrix builder.
     * @param k wavenumber (complex)
     * @param c_o refraction index of outer domain (should be 1)
     * @param c_i refraction indef of inner domain (must not be smaller than c_o)
     * @param T complex matrix to which the solutions operator matrix will be stored
     * @param T_der complex matrix to which the 1st derivative of the solutions operator matrix will be stored
     * @param T_der2 complex matrix to which the 2nd derivative of the solutions operator matrix will be stored
     */
    void gen_sol_op_2nd_der(const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2);
    /**
     * Compute approximation of solutions operator and its 1st and 2nd derivatives
     * for second-kind direct BIEs of Helmholtz transmission problem using Galerkin BEM.
     * This routine will use the given Galerkin matrix builder.
     * @param builder Galerkin matrix builder
     * @param k wavenumber (complex)
     * @param c_o refraction index of outer domain (should be 1)
     * @param c_i refraction indef of inner domain (must not be smaller than c_o)
     * @param T complex matrix to which the solutions operator matrix will be stored
     * @param T_der complex matrix to which the 1st derivative of the solutions operator matrix will be stored
     * @param T_der2 complex matrix to which the 2nd derivative of the solutions operator matrix will be stored
     */
    void gen_sol_op_2nd_der(GalerkinMatrixBuilder &builder, const std::complex<double> &k, double c_o, double c_i,
                            Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2);
    /**
     * Return reference to the mass matrix.
     */
    const Eigen::MatrixXd &mass_matrix() const { return M; };
};

#endif //GEN_SOL_OPHPP
