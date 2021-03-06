/**
 * \file double_layer_der2.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        second derivative in \p k of the Double Layer BIO for the
 *        Helmholtz kernel, using
 *        common and composite Gauss-Legendre quadrature rules.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#ifndef DOUBLELAYERDER2HPP
#define DOUBLELAYERDER2HPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"

/**
 * \namespace double_layer_helmholtz_der2
 * \brief This namespace contains all the functions for evaluating the second
 * derivative in \p k of the Double Layer Galerkin Matrix using
 * common and composite Gauss-Legendre quadrature rules.
 */
    namespace double_layer_helmholtz_der2 {
/**
 * This function is used to evaluate the Interaction Matrix for a pair of
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the second
 * derivative in \p k of the Double Layer BIO for the Helmholtz Kernel,
 * using the given trial and test spaces.
 * It implements the case where the panels \f$\Pi\f$ and \f$\Pi\f$' are adjacent.
 * This function calculates a matrix entry by using a composite Gauss-Legnedre
 * quadrature rule with which a tensor product quadrature rule is generated.
 *
 * @param pi parametrization for the first panel \f$\Pi\f$
 * @param pi_p parametrization for the second panel \f$\Pi\f$'
 * @param trial_space the trial space for evaluating the matrix
 * @param test_space the test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied
 * @param k wavenumber
 * @param c refraction index
 * @return the matrix K'' for the second derivative in \p k of the
 * Double Layer BIO bilinear form
 */
        Eigen::MatrixXcd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                                 const AbstractParametrizedCurve &pi_p,
                                                 const AbstractBEMSpace &trial_space,
                                                 const AbstractBEMSpace &test_space,
                                                 const QuadRule &GaussQR,
                                                 const std::complex<double> k,
                                                 const double c);

/**
 * This function is used to evaluate the Interaction Matrix for a pair of
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the
 * second derivative in \p k of the Double Layer BIO for the Helmholtz kernel,
 * using the given trial and test spaces.
 * It implements the case where the panels \f$\Pi\f$ and \f$\Pi\f$' are coinciding.
 * This function calculates a matrix entry by using a composite Gauss-Legnedre
 * quadrature rule with which a tensor product quadrature rule is generated.
 *
 * @param pi parametrization for the first panel \f$\Pi\f$
 * @param pi_p parametrization for the second panel \f$\Pi\f$'
 * @param trial_space the trial space for evaluating the matrix
 * @param test_space the test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied
 * @param k wavenumber
 * @param c refraction index
 * @return the matrix K'' for the second derivative in \p k of the
 * Double Layer BIO bilinear form
 */
        Eigen::MatrixXcd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                                   const AbstractParametrizedCurve &pi_p,
                                                   const AbstractBEMSpace &trial_space,
                                                   const AbstractBEMSpace &test_space,
                                                   const QuadRule &GaussQR,
                                                   const std::complex<double> k,
                                                   const double c);

/**
 * This function is used to evaluate the Interaction Matrix for a pair of
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the
 * second derivative in \p k of the Double Layer BIO for the Helmholtz kernel,
 * using the given trial and test spaces.
 * It implements the case where the panels \f$\Pi\f$ and \f$\Pi\f$' are
 * completely disjoint. This function calculates a matrix entry by using
 * Gauss Legendre quadrature rule.
 *
 * @param pi parametrization for the first panel \f$\Pi\f$
 * @param pi_p parametrization for the second panel \f$\Pi\f$'
 * @param trial_space the trial space for evaluating the matrix
 * @param test_space the test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied
 * @param k wavenumber
 * @param c refraction index
 * @return the matrix K'' for the second derivative in \p k of the
 * Double Layer BIO bilinear form
 */
        Eigen::MatrixXcd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                                const AbstractParametrizedCurve &pi_p,
                                                const AbstractBEMSpace &trial_space,
                                                const AbstractBEMSpace &test_space,
                                                const QuadRule &GaussQR,
                                                const std::complex<double> k,
                                                const double c);

/**
 * This function is used to evaluate the Interaction Matrix
 * for the pair of panels \f$\Pi\f$ and \f$\Pi\f$' for the
 * bilinear form induced by the second derivative in \p k of
 * the Double Layer BIO for the Helmholtz kernel, given by the formula:
 * \f{eqnarray*}{
 * I_{ij} = \frac{ic}{4} \int_{0}^{1} &\int_{0}^{1}&
 * &\Big(& k\sqrt{c}H_1''^{(1)}(k\sqrt{c}\|\gamma_{\Pi}(s)-\gamma_{\Pi'}(t)\|)
 * \|\gamma_{\Pi}(s)-\gamma_{\Pi'}(t)\| \\
 * &&&&+ 2H_1'^{(1)}(k\sqrt{c}\|\gamma_{\Pi}(s)-\gamma_{\Pi'}(t)\|) \Big) \\
 * &&&\cdot&(\gamma_{\Pi}(s)-\gamma_{\Pi'}(t))
 * \cdot \textbf{n}(\gamma_{\Pi'}(t)) 
 * \hat{b}^{j}(t) \hat{\beta}^{i}(s) \|\dot{\gamma}_{\Pi}(s)\|
 * \|\dot{\gamma}_{\Pi'}(t)\| dt ds 
 * \f} 
 * is computed.
 * For \f$ k \rightarrow 0 \f$ the limit
 * \f$I_{ij}\f$ = 0
 * is computed since in the limit the double layer kernel does not depend on \f$k\f$.
 * \f$I\f$, the interaction matrix is of size \f$Q_{test}\times Q_{trial}\f$
 * where \f$Q_{test}\f$ is the number of reference shape functions for the test
 * BEM space and \f$Q_{trial}\f$ is the number of reference shape functions in
 * the trial BEM space. The computation of the entries are based on cases and
 * delegated to these functions accordingly:
 *
 * ComputeIntegralGeneral()
 *
 * ComputeIntegralAdjacent()
 *
 * ComputeIntegralCoinciding()
 *
 * @param pi parametrization for the first panel \f$\Pi\f$
 * @param pi_p parametrization for the second panel \f$\Pi\f$'
 * @param trial_space the trial space for evaluating the matrix
 * @param test_space the test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the
 * Gauss-Legendre quadrature rule to be applied.
 * @param CGaussQR QuadRule object containing the composite
 * Gauss-Legendre quadrature rule to be applied.
 * @param k wavenumber
 * @param c refraction index
 * @return An Eigen::MatrixXd type Interaction Matrix
 * (\f$Q_{test}\times Q_{trial}\f$)
 */
        Eigen::MatrixXcd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                           const AbstractParametrizedCurve &pi_p,
                                           const AbstractBEMSpace &trial_space,
                                           const AbstractBEMSpace &test_space,
                                           const QuadRule &GaussQR,
                                           const QuadRule &CGaussQR,
                                           const std::complex<double> k,
                                           const double c);

/**
 * This function is used to evaluate the full Galerkin matrix based on the
 * second derivative in \p k of the bilinear form for Double Layer BIO for
 * the Helmholtz kernel. It uses the trial and test spaces
 * and the parametrized mesh object, specified in the inputs to the function.
 * It evaluates the matrix by panel oriented assembly by
 * first evaluating the interaction matrix for all possible pairs of panels and
 * then using the local to global map of BEM spaces to fill the matrix entries.
 *
 * @param mesh ParametrizedMesh object containing all the panels in the form
 *             of small parametrized curves
 * @param trial_space the trial space for evaluating the matrix
 * @param test_space the test space for evaluating the matrix
 * @param N the order for gauss/log-weighted quadrature
 * @param k wavenumber
 * @param c refraction index
 * @return an Eigen::MatrixXd type Galerkin Matrix for the given mesh and space
 */
        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &trial_space,
                                        const AbstractBEMSpace &test_space,
                                        const unsigned int &N,
                                        const std::complex<double> k,
                                        const double c);

    }

#endif // DOUBLELAYERHPP
