/**
 * \file galerkin_matrix_builder.hpp
 * \brief This file declares the object to evaluate the entries of
 *        Galerkin matrices by using common and composite
 *        Gauss-Legendre quadrature rules.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 * (c) 2023 Luka Marohnić
 */

#ifndef GALERKIN_ALLHPP
#define GALERKIN_ALLHPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"

typedef std::complex<double> complex_t;

class BuilderData
{
    // panel vector
    PanelVector panels;
    // dimensions
    size_t numpanels, Qtest, Qtrial, dim_test, dim_trial;
public:
    // trial and test spaces
    const ParametrizedMesh &mesh;
    const AbstractBEMSpace &test_space;
    const AbstractBEMSpace &trial_space;
    // quadrature rules
    QuadRule GaussQR;
    QuadRule CGaussQR;
    /**
     * Construct the Galerkin builder data object from the given parameters.
     * @param mesh panel mesh
     * @param test_space_in the test space
     * @param trial_space_in the trial space
     * @param order the order of Gaussian quadrature
     */
    BuilderData(const ParametrizedMesh &mesh_in,
                const AbstractBEMSpace &test_space_in,
                const AbstractBEMSpace &trial_space_in,
                unsigned order);
    /**
     * Return the dimension of the test space.
     * @return unsigned int, the dimension
     */
    size_t getTestSpaceDimension() const { return dim_test; }
    /**
     * Return the dimension of the trial space.
     * @return unsigned int, the dimension
     */
    size_t getTrialSpaceDimension() const { return dim_trial; }
    /**
     * Return TRUE iff test and trial spaces are equal (i.e. the same space).
     * @return boolean true if spaces are equal, false otherwise
     */
    bool testTrialSpacesAreEqual() const { return &test_space == &trial_space; }
    /**
     * Return TRUE iff i-th and j-th panels are adjacent.
     * @return boolean true if panels are adjacent, false otherwise
     */
    bool is_adjacent(size_t i, size_t j, bool &swap) const;
    // data
    Eigen::ArrayXXd m_tc, m_ta, m_tg, m_sc, m_sa, m_sg, m_wc, m_wa, m_wg;
    Eigen::ArrayXXd m_double_layer_coinciding_fg, m_double_layer_adjacent_fg, m_double_layer_adjacent_fg_swap, m_double_layer_general_fg,
                    m_double_layer_coinciding_fg_t, m_double_layer_adjacent_fg_t, m_double_layer_adjacent_fg_swap_t,
                    m_double_layer_general_fg_t;
    Eigen::ArrayXXd m_single_layer_coinciding_fg, m_single_layer_adjacent_fg, m_single_layer_adjacent_fg_swap, m_single_layer_general_fg;
    Eigen::ArrayXXd m_hypersingular_coinciding_fg, m_hypersingular_adjacent_fg, m_hypersingular_adjacent_fg_swap, m_hypersingular_general_fg,
                    m_hypersingular_coinciding_fg_arc, m_hypersingular_adjacent_fg_arc, m_hypersingular_adjacent_fg_arc_swap,
                    m_hypersingular_general_fg_arc;
    Eigen::ArrayXXcd Derivative_01_sc, Derivative_01_tc, Derivative_01_sa, Derivative_01_ta, Derivative_01_sg, Derivative_01_tg;
    Eigen::ArrayXXcd Derivative_01_swapped_sa, Derivative_01_swapped_ta;
    Eigen::ArrayXXd Derivative_01_sc_n, Derivative_01_tc_n, Derivative_01_sa_n, Derivative_01_ta_n, Derivative_01_sg_n, Derivative_01_tg_n;
    Eigen::ArrayXXd Derivative_01_swapped_sa_n, Derivative_01_swapped_ta_n;
    Eigen::ArrayXXcd op_sc, op_tc, op_sa, op_ta, op_sg, op_tg, swapped_op_sa, swapped_op_ta;
    // methods
    size_t getQtest() const { return Qtest; }
    size_t getQtrial() const { return Qtrial; }
    size_t getGaussQROrder() const { return GaussQR.n; }
    size_t getCGaussQROrder() const { return CGaussQR.n; }
};

class GalerkinMatrixBuilder
{
    // builder data object
    const BuilderData &data;
    // wavenumber and refraction index
    complex_t k, ksqrtc, kkc;
    double c, sqrtc, ksqrtca;
    bool k_real_positive;
    // timing info
    unsigned interaction_matrix_assembly_time;
    unsigned hankel_computation_time;
    // workspace
    Eigen::ArrayXXd m_zero, m_zero_s;
    Eigen::ArrayXXcd m_h0[2], m_h1[2], m_v[2], m_tangent[2], m_tangent_p[2], masked1[2], masked2[2], masked3[2];
    Eigen::ArrayXXcd m_h0_s, m_h1_s, m_v_s, m_tangent_s, m_tangent_p_s, masked1_s, masked2_s, masked3_s;
    Eigen::ArrayXXd m_v_norm[2], m_v_norm2[2], m_tangent_norm[2], m_tangent_p_norm[2];
    Eigen::ArrayXXd m_v_norm_s, m_v_norm2_s, m_tangent_norm_s, m_tangent_p_norm_s;
    // interaction matrices
    Eigen::ArrayXXcd double_layer_interaction_matrix, hypersingular_interaction_matrix, single_layer_interaction_matrix;
    Eigen::ArrayXXcd double_layer_der_interaction_matrix, hypersingular_der_interaction_matrix, single_layer_der_interaction_matrix;
    Eigen::ArrayXXcd double_layer_der2_interaction_matrix, hypersingular_der2_interaction_matrix, single_layer_der2_interaction_matrix;
    // assembled matrices
    Eigen::MatrixXcd double_layer_matrix, hypersingular_matrix, single_layer_matrix;
    Eigen::MatrixXcd double_layer_der_matrix, hypersingular_der_matrix, single_layer_der_matrix;
    Eigen::MatrixXcd double_layer_der2_matrix, hypersingular_der2_matrix, single_layer_der2_matrix;
    // routines for computing shared data
    void compute_coinciding(size_t i) throw();
    void compute_adjacent(size_t i, size_t j, bool swap) throw();
    void compute_general(size_t i, size_t j) throw();
    // routines for interaction matrix construction
    void double_layer_coinciding(int der) throw();
    void double_layer_adjacent(bool swap, int der, bool transp) throw();
    void double_layer_general(int der, bool transp) throw();
    void hypersingular_coinciding(int der) throw();
    void hypersingular_adjacent(bool swap, int der) throw();
    void hypersingular_general(int der) throw();
    void single_layer_coinciding(int der) throw();
    void single_layer_adjacent(bool swap, int der) throw();
    void single_layer_general(int der) throw();
    void all_coinciding(int der) throw();
    void all_adjacent(bool swap, int der) throw();
    void all_general(int der) throw();
    void initialize_parameters(const std::complex<double> &k_in, double c_in);

public:
    /**
     * Construct the Galerkin matrix builder from the given parameters.
     * @param builder_data the Galerkin builder data object
     */
    GalerkinMatrixBuilder(const BuilderData &builder_data);
    /**
     * Get the DER-th derivative of the previously assembled K submatrix.
     * @param der the order of derivative (default 0)
     * @return the reference to the assembled matrix
     */
    const Eigen::MatrixXcd &getDoubleLayer(int der = 0) const;
    /**
     * Get the DER-th derivative of the previously assembled W submatrix.
     * @param der the order of derivative (default 0)
     * @return the reference to the assembled matrix
     */
    const Eigen::MatrixXcd &getHypersingular(int der = 0) const;
    /**
     * Get the DER-th derivative of the previously assembled V submatrix.
     * @param der the order of derivative (default 0)
     * @return the reference to the assembled matrix
     */
    const Eigen::MatrixXcd &getSingleLayer(int der = 0) const;
    /**
     * Assemble the K submatrix and optionally its first DER derivatives
     * for the given wavenumber and refraction index.
     * @param k_in wavenumber (complex)
     * @param c_in refraction index (real, at least 1 otherwise an error is thrown)
     * @param der the order of derivative (default 0, derivatives are not computed)
     */
    void assembleDoubleLayer(const std::complex<double> &k_in, double c_in, int der = 0);
    /**
     * Assemble the W submatrix and optionally its first DER derivatives
     * for the given wavenumber and refraction index.
     * @param k_in wavenumber (complex)
     * @param c_in refraction index (real, at least 1 otherwise an error is thrown)
     * @param der the order of derivative (default 0, derivatives are not computed)
     */
    void assembleHypersingular(const std::complex<double> &k_in, double c_in, int der = 0);
    /**
     * Assemble the V submatrix and optionally its first DER derivatives
     * for the given wavenumber and refraction index.
     * @param k_in wavenumber (complex)
     * @param c_in refraction index (real, at least 1 otherwise an error is thrown)
     * @param der the order of derivative (default 0, derivatives are not computed)
     */
    void assembleSingleLayer(const std::complex<double> &k_in, double c_in, int der = 0);
    /**
     * Assemble K, W and V submatrices and optionally their first DER derivatives
     * for the given wavenumber and refraction index. An error will be thrown
     * if the test and trial spaces are not equal.
     * @param k_in wavenumber (complex)
     * @param c_in refraction index (real, at least 1 otherwise an error is thrown)
     * @param der the order of derivative (default 0, derivatives are not computed)
     */
    void assembleAll(const std::complex<double> &k_in, double c_in, int der = 0);
    // get timing info for profiling
    unsigned getInteractionMatrixAssemblyTime();
    unsigned getHankelComputationTime();
};

#endif // GALERKIN_ALLHPP
