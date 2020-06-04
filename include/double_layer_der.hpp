/**
 * \file double_layer.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Double Layer BIO, using the transformations given in
 *        \f$\ref{ss:quadapprox}\f$ in the Lecture Notes for Advanced Numerical
 *        Methods for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef DOUBLELAYERDERHPP
#define DOUBLELAYERDERHPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"

/**
 * This namespace contains all the functions for evaluating the Double Layer
 * Galerkin Matrix using quadrature and panel oriented assembly
 */
    namespace double_layer_helmholtz_der {
/**
 * This function is used to evaluate the Interaction Matrix for a pair of
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the Double
 * Layer BIO, using the given trial and test spaces. It implements the case
 * where the panels \f$\Pi\f$ and \f$\Pi\f$' are adjacent. This function
 * calculates a matrix entry by using a local arclength parametrization given in
 * \f$\eqref{eq:ap}\f$ and the transformations mentioned in
 * \f$\ref{par:Kadpan}\f$
 *
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param trial_space The trial space for evaluating the matrix.
 * @param test_space The test space for evaluating the matrix.
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied.
 * @return The matrix K for Double Layer BIO bilinear form.
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
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the Double
 * Layer BIO, using the given trial and test spaces. It implements the case
 * where the panels \f$\Pi\f$ and \f$\Pi\f$' are coinciding. This function
 * calculates a matrix entry by using the transformations mentioned in
 * \f$\ref{par:Kidpan}\f$
 *
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param trial_space The trial space for evaluating the matrix.
 * @param test_space The test space for evaluating the matrix.
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied.
 * @return The matrix K for Double Layer BIO bilinear form.
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
 * panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear form induced by the Double
 * Layer BIO, using the given trial and test spaces. It implements the case
 * where the panels \f$\Pi\f$ and \f$\Pi\f$' are completely disjoint. This
 * function calculates a matrix entry by using Gauss Legendre quadrature rule.
 *
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param trial_space The trial space for evaluating the matrix.
 * @param test_space The test space for evaluating the matrix.
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied.
 * @return The matrix K for Double Layer BIO bilinear form.
 */
        Eigen::MatrixXcd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                                const AbstractParametrizedCurve &pi_p,
                                                const AbstractBEMSpace &trial_space,
                                                const AbstractBEMSpace &test_space,
                                                const QuadRule &GaussQR,
                                                const std::complex<double> k,
                                                const double c);

/**
 * This function is used to evaluate the Interaction Matrix defined in
 * \f$\eqref{eq:Al}\f$ for the pair of panels \f$\Pi\f$ and \f$\Pi\f$' for the
 * bilinear form induced by the Double Layer BIO; given by the formula :
 * \f$I_{ij}\f$ = \f$-\frac{1}{2\pi} \int_{-1}^{1} \int_{-1}^{1}
 * \frac{(\gamma_{\Pi}(s)-\gamma_{\Pi'}(t))}
 * {\|\gamma_{\Pi}(s)-\gamma_{\Pi'}(t)\|^2}.\textbf{n}(\gamma_{\Pi'}(t))
 * \hat{b}^{j}(t) \hat{\beta}^{i}(s) \|\dot{\gamma}_{\Pi}(s)\|
 * \|\dot{\gamma}_{\Pi'}(t)\| dt ds \f$ where \f$\hat{b}^{j}\f$
 *  & \f$\hat{\beta}^{i}\f$ are reference shape functions associated with the
 * trial space \f$S_{p}^{0}\f$ and test space \f$S_{p}^{-1}\f$ respectively.
 * \f$I\f$, the interaction matrix is of size \f$Q_{test}\f$X\f$Q_{trial}\f$
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
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param trial_space The trial space for evaluating the matrix.
 * @param test_space The test space for evaluating the matrix.
 * @param GaussQR QuadRule object containing the Gaussian quadrature to be
 * applied.
 * @return An Eigen::MatrixXd type Interaction Matrix
 * (\f$Q_{test}\f$X\f$Q_{trial}\f$)
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
 * Bilinear form for Double Layer BIO. It uses the trial and test spaces
 * and the parametrized mesh object, specified in the inputs to the function.
 * It evaluates the matrix by panel oriented assembly (\f$\ref{pc:ass}\f$) by
 * first evaluating the interaction matrix for all possible pairs of panels and
 * then using the local to global map of BEM spaces to fill the matrix entries.
 *
 * @param mesh ParametrizedMesh object containing all the panels in the form
 *             of small parametrized curves
 * @param trial_space The trial space for evaluating the matrix.
 * @param test_space The test space for evaluating the matrix.
 * @param N The order for gauss/log-weighted quadrature.
 * @return An Eigen::MatrixXd type Galerkin Matrix for the given mesh and space
 */
        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &trial_space,
                                        const AbstractBEMSpace &test_space,
                                        const unsigned int &N,
                                        const std::complex<double> k,
                                        const double c);

    }

#endif // DOUBLELAYERHPP
