/**
 * \file hypersingular.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Hypersingular BIO for the Helmholtz kernel,
 *        using common and composite Gauss-Legendre quadrature rules.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#ifndef HYPERSINGULARHPP
#define HYPERSINGULARHPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"

/**
 * \namespace hypersingular_helmholtz 
 * \brief This namespace contains all the functions for evaluating the Hypersingular
 * Galerkin Matrix using common and composite Gauss-Legendre quadrature rules.
 */
namespace hypersingular_helmholtz {

/**
 * This function is used to evaluate the Interaction Matrix for a pair of
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the
 * Hypersingular BIO for the Helmholtz kernel, using the given trial and test spaces.
 * It implements the case where the panels \f$\Pi\f$ and \f$\Pi\f$' are adjacent.
 * This function calculates a matrix entry by using a composite Gauss-Legnedre
 * quadrature rule with which a tensor product quadrature rule is generated.
 *
 * @param pi parametrization for the first panel \f$\Pi\f$
 * @param pi_p parametrization for the second panel \f$\Pi\f$'
 * @param space the trial & test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied
 * @param k wavenumber
 * @param c refraction index
 * @return the matrix W for Helmholtz kernel Hypersingular BIO bilinear form
 */
    void ComputeIntegralAdjacent(Eigen::MatrixXcd &interaction_matrix,
                                 const AbstractParametrizedCurve &pi,
                                 const AbstractParametrizedCurve &pi_p,
                                 const AbstractBEMSpace &space,
                                 const QuadRule &GaussQR,
                                 const std::complex<double> &k,
                                 const double c_i, const double c_o,
                                 gq_workspace_t &ws);

/**
 * This function is used to evaluate the Interaction Matrix for a pair of
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the
 * Hypersingular BIO for the Helmholtz kernel, using the given trial and test spaces.
 * It implements the case where the panels \f$\Pi\f$ and \f$\Pi\f$' are coinciding.
 * This function calculates a matrix entry by using a composite Gauss-Legnedre
 * quadrature rule with which a tensor product quadrature rule is generated.
 *
 * @param pi parametrization for the first panel \f$\Pi\f$
 * @param pi_p parametrization for the second panel \f$\Pi\f$'
 * @param space the trial & test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied
 * @param k wavenumber
 * @param c refraction index
 * @return the matrix W for Helmholtz kernel Hypersingular BIO bilinear form
 */
    void ComputeIntegralCoinciding(Eigen::MatrixXcd &interaction_matrix,
                                   const AbstractParametrizedCurve &pi,
                                   const AbstractParametrizedCurve &pi_p,
                                   const AbstractBEMSpace &space,
                                   const QuadRule &GaussQR,
                                   const std::complex<double> &k,
                                   const double c_i, const double c_o,
                                   gq_workspace_t &ws);

/**
 * This function is used to evaluate the Interaction Matrix for a pair of
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the
 * Hypersingular BIO for the Helmholtz kernel, using the given trial and test spaces. It implements the case
 * where the panels \f$\Pi\f$ and \f$\Pi\f$' are completely disjoint. This
 * function calculates a matrix entry by using Gauss-Legendre quadrature rule.
 *
 * @param pi parametrization for the first panel \f$\Pi\f$
 * @param pi_p parametrization for the second panel \f$\Pi\f$'
 * @param space the trial & test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied
 * @param k wavenumber
 * @param c refraction index
 * @return the matrix W for Helmholtz kernel Hypersingular BIO bilinear form
 */
    void ComputeIntegralGeneral(Eigen::MatrixXcd &interaction_matrix,
                                const AbstractParametrizedCurve &pi,
                                const AbstractParametrizedCurve &pi_p,
                                const AbstractBEMSpace &space,
                                const QuadRule &GaussQR,
                                const std::complex<double> &k,
                                const double c_i, const double c_o,
                                gq_workspace_t &ws);

/**
 * This function is used to evaluate the Interaction Matrix
 * for the pair of panels \f$\Pi\f$ and \f$\Pi\f$' for the
 * bilinear form induced by the Helmholtz kernel Hypersingular BIO,
 * given the formula:
 * \f{eqnarray*}{ 
 * I_{ij} = \frac{i}{4} \int_{0}^{1} \int_{0}^{1}
 * &H_0^{(1)}&(k\sqrt{c}\|\gamma_{\Pi}(s)-\gamma_{\Pi'}(t)\|) \\
 * &\cdot& \Big[\hat{B}^{j}(t) \hat{B}^{i}(s) 
 * - ck^2 \hat{b}^{j}(t) \hat{b}^{i}(s) \|\dot{\gamma}_{\Pi}(s)\| \|\dot{\gamma}_{\Pi'}(t)\|
 * \mathbf{n}(\gamma_{\Pi}(s)) \cdot \mathbf{n}(\gamma_{\Pi'}(s))\Big]
 * dt ds \f}
 * where \f$\hat{B}^{j}\f$ and \f$\hat{B}^{j}\f$
 * refer to the arc-length derivatives of \f$\hat{b}^{j}\f$ and \f$\hat{b}^{j}\f$ which 
 * are reference shape functions associated with the
 * trial space and test space \f$S_{p}^{1}\f$.
 * For \f$ k \rightarrow 0 \f$ the limit
 * \f{eqnarray*}{
 * I_{ij} = -\frac{1}{2\pi} \int_{0}^{1} \int_{0}^{1}
 * &\ln&(\|\gamma_{\Pi}(s)-\gamma_{\Pi'}(t)\|) \\
 * &\cdot& \Big[\hat{B}^{j}(t) \hat{B}^{i}(s) 
 * - (ck^2)\hat{b}^{j}(t) \hat{b}^{i}(s) \|\dot{\gamma}_{\Pi}(s)\| \|\dot{\gamma}_{\Pi'}(t)\|
 * \mathbf{n}(\gamma_{\Pi}(s)) \cdot \mathbf{n}(\gamma_{\Pi'}(s))\Big]
 * dt ds 
 * \f}
 * is computed.
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
 * @param space the trial & test space for evaluating the matrix
 * @param GaussQR QuadRule object containing the
 * Gauss-Legendre quadrature rule to be applied.
 * @param CGaussQR QuadRule object containing the composite
 * Gauss-Legendre quadrature rule to be applied.
 * @param k wavenumber
 * @param c refraction index
 * @return an Eigen::MatrixXd type Interaction Matrix
 * (\f$Q_{test}\times Q_{trial}\f$)
 */
    void InteractionMatrix(Eigen::MatrixXcd &interaction_matrix,
                           const AbstractParametrizedCurve &pi,
                           const AbstractParametrizedCurve &pi_p,
                           const AbstractBEMSpace &space,
                           const QuadRule &GaussQR,
                           const QuadRule &CGaussQR,
                           const std::complex<double> &k,
                           const double c_i, const double c_o,
                           gq_workspace_t &ws);

/**
 * This function is used to evaluate the full Galerkin matrix based on the
 * Bilinear form for the Helmholtz kernel Hypersingular BIO. It uses the trial
 * and test spaces and the parametrized mesh object, specified in the inputs
 * to the function. It evaluates the matrix by panel oriented assembly by
 * first evaluating the interaction matrix for all possible pairs of panels and
 * then using the local to global map of BEM spaces to fill the matrix entries.
 *
 * @param mesh ParametrizedMesh object containing all the panels in the form
 *             of small parametrized curves
 * @param space the trial & test space for evaluating the matrix
 * @param N the order for common/composite Gauss-Legendre quadrature rule
 * @param k wavenumber
 * @param c refraction index
 * @return an Eigen::MatrixXd type Galerkin Matrix for the given mesh and space
 */
    Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                    const AbstractBEMSpace &space,
                                    const unsigned int &N,
                                    const std::complex<double> &k,
                                    const double c_i,
                                    const double c_o);


} // namespace hypersingular_helmholtz

#endif // HYPERSINGULARHPP
