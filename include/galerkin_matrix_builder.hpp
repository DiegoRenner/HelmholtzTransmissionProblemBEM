/**
 * \file galerkin_matrix_builder.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices by using common and composite
 *        Gauss-Legendre quadrature rules.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#ifndef GALERKIN_ALLHPP
#define GALERKIN_ALLHPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"

typedef std::complex<double> complex_t;

class GalerkinMatrixBuilder
{
    // trial and test spaces
    const AbstractBEMSpace &test_space;
    const AbstractBEMSpace &trial_space;
    // panel vector
    const PanelVector &panels;
    // quadrature rules
    const QuadRule &GaussQR;
    const QuadRule &CGaussQR;
    // wavenumber and refraction index
    complex_t k, ksqrtc, kkc;
    double c, sqrtc, ksqrtca;
    // dimensions
    size_t numpanels, Qtest, Qtrial, dim_test, dim_trial;
    // workspace
    Eigen::MatrixXcd m_h0, m_h1, m_v, m_tangent, m_tangent_p;
    Eigen::MatrixXd m_v_norm, m_tangent_norm, m_tangent_p_norm;
    std::vector<size_t> ind2N;
    // interaction matrices
    Eigen::MatrixXcd double_layer_interaction_matrix, hypersingular_interaction_matrix, single_layer_interaction_matrix;
    Eigen::MatrixXcd double_layer_der_interaction_matrix, hypersingular_der_interaction_matrix, single_layer_der_interaction_matrix;
    Eigen::MatrixXcd double_layer_der2_interaction_matrix, hypersingular_der2_interaction_matrix, single_layer_der2_interaction_matrix;
    // assembled matrices
    Eigen::MatrixXcd double_layer_matrix, hypersingular_matrix, single_layer_matrix;
    Eigen::MatrixXcd double_layer_der_matrix, hypersingular_der_matrix, single_layer_der_matrix;
    Eigen::MatrixXcd double_layer_der2_matrix, hypersingular_der2_matrix, single_layer_der2_matrix;
    // routines for computing shared data
    void compute_coinciding(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p) throw();
    void compute_adjacent(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p, bool swap) throw();
    void compute_general(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p) throw();
    // routines for interaction matrix construction
    void double_layer_coinciding(int der) throw();
    void double_layer_adjacent(bool swap, int der) throw();
    void double_layer_general(int der) throw();
    void hypersingular_coinciding(int der) throw();
    void hypersingular_adjacent(bool swap, int der) throw();
    void hypersingular_general(int der) throw();
    void single_layer_coinciding(int der) throw();
    void single_layer_adjacent(bool swap, int der) throw();
    void single_layer_general(int der) throw();
    void initialize_parameters(const std::complex<double> &k_in, double c_in);
    bool is_adjacent(const AbstractParametrizedCurve &p1, const AbstractParametrizedCurve &p2, bool &swap) const;

public:
    GalerkinMatrixBuilder(const ParametrizedMesh &mesh,
                          const AbstractBEMSpace &test_space_in,
                          const AbstractBEMSpace &trial_space_in,
                          const QuadRule &GaussQR_in,
                          const QuadRule &CGaussQR_in);
    size_t getTestSpaceDimension() const { return dim_test; }
    size_t getTrialSpaceDimension() const { return dim_trial; }
    bool testTrialSpacesAreEqual() const { return &test_space == &trial_space; }
    const Eigen::MatrixXcd &getDoubleLayer(int der = 0) const;
    const Eigen::MatrixXcd &getHypersingular(int der = 0) const;
    const Eigen::MatrixXcd &getSingleLayer(int der = 0) const;
    // assembly routines
    void assembleDoubleLayer(const std::complex<double> &k_in, double c_in, int der = 0);
    void assembleHypersingular(const std::complex<double> &k_in, double c_in, int der = 0);
    void assembleSingleLayer(const std::complex<double> &k_in, double c_in, int der = 0);
    void assembleAll(const std::complex<double> &k_in, double c_in, int der = 0);
};


#endif // GALERKIN_ALLHPP

