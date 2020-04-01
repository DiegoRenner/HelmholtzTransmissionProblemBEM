/**
 * \file single_layer.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Single Layer BIO, using the transformations given in
 *        \f$\ref{ss:quadapprox}\f$ in the Lecture Notes for Advanced Numerical
 *        Methods for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef SINGLELAYERCORRECTIONHPP
#define SINGLELAYECORRECTIONRHPP

#include "abstract_bem_space.hpp"
#include "logweight_quadrature.hpp"

namespace parametricbem2d {
/**
 * This namespace contains all the functions for evaluating the Single Layer
 * Galerkin Matrix using quadrature and panel oriented assembly
 */
namespace single_layer_correction {

/**
 * This function is used to evaluate the Interaction Matrix for the pair
 * of panels \f$\Pi\f$ and \f$\Pi\f$', for the bilinear form induced by
 * the Single Layer BIO. It implements the case where the panels are adjacent.
 * The matrix entries are calculated by using a local arclength parametrization
 * and the transformations mentioned in \f$\ref{par:Vadj}\f$
 *
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param space The BEM space to be used for calculations
 * @param GaussQR QuadRule object containing the Gaussian Quadrature to be
 * applied.
 * @return An Eigen::MatrixXd type Interaction Matrix (QXQ)
 *         where Q is number of local shape functions in BEM space
 */
Eigen::MatrixXd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                        const AbstractParametrizedCurve &pi_p,
                                        const AbstractBEMSpace &space,
                                        const QuadRule &GaussQR,
                                        const double k);

/**
 * This function is used to evaluate the Interaction Matrix for the pair
 * of panels \f$\Pi\f$ and \f$\Pi\f$', for the bilinear form induced by
 * the Single Layer BIO. It implements the case where the panels are coinciding.
 * The matrix entries are calculated by using transformations mentioned in
 * \f$\ref{par:Vidp}\f$
 *
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param space The BEM space to be used for calculations
 * @param GaussQR QuadRule object containing the Gaussian Quadrature to be
 * applied.
 * @return An Eigen::MatrixXd type Interaction Matrix (QXQ)
 *         where Q is number of local shape functions in BEM space
 */
Eigen::MatrixXd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                          const AbstractParametrizedCurve &pi_p,
                                          const AbstractBEMSpace &space,
                                          const QuadRule &GaussQR,
                                          const double k);

/**
 * This function is used to evaluate the Interaction Matrix for the pair
 * of panels \f$\Pi\f$ and \f$\Pi\f$', for the bilinear form induced by
 * the Single Layer BIO. It implements the case where the panels are completely
 * disjoint. The matrix entries are calculated by Gauss Legendre quadrature.
 *
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param space The BEM space to be used for calculations
 * @param GaussQR QuadRule object containing the Gaussian Quadrature to be
 * applied.
 * @return An Eigen::MatrixXd type Interaction Matrix (QXQ)
 *         where Q is number of local shape functions in BEM space
 */
Eigen::MatrixXd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                       const AbstractParametrizedCurve &pi_p,
                                       const AbstractBEMSpace &space,
                                       const QuadRule &GaussQR,
                                       const double k);
/**
 * This function is used to evaluate the Interaction Matrix defined in
 * \f$\eqref{eq:Al}\f$ for the pair of panels \f$\Pi\f$ and \f$\Pi\f$', for the
 * bilinear form induced by the Single Layer BIO given by the formula :
 * \f$I_{ij}\f$ = \f$-\frac{1}{2\pi} \int_{-1}^{1} \int_{-1}^{1}
 * \log{(\|\gamma_{\Pi}(s)-\gamma_{\Pi'}(t)\|)} \hat{b}^{j}(t) \hat{b}^{i}(s)
 * \|\dot{\gamma}_{\Pi}(s)\| \|\dot{\gamma}_{\Pi'}(t)\| dt ds \f$ where
 * \f$\hat{b}^{j}\f$ & \f$\hat{b}^{i}\f$ are local shape functions associated
 * the trial and test BEM space \f$S_{p}^{-1}\f$. The interaction matrix \f$I\f$
 * , is of size QXQ where Q is the number of local shape functions in the BEM
 * space. The computation of the entries are based on cases and delegated to
 * these functions accordingly:
 *
 * ComputeIntegralGeneral()
 *
 * ComputeIntegralAdjacent()
 *
 * ComputeIntegralCoinciding()
 *
 * @param pi Parametrization for the first panel \f$\Pi\f$.
 * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
 * @param space The BEM space to be used for calculations
 * @param GaussQR QuadRule object containing the Gaussian Quadrature to be
 * applied.
 * @return An Eigen::MatrixXd type Interaction Matrix (QXQ)
 *         where Q is number of local shape functions in BEM space
 */
Eigen::MatrixXd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                  const AbstractParametrizedCurve &pi_p,
                                  const AbstractBEMSpace &space,
                                  const QuadRule &GaussQR,
                                  const double k);

/**
 * This function is used to evaluate the full Galerkin matrix based on the
 * Bilinear form for Single Layer BIO. It uses the trial and test spaces and
 * Parametrized mesh specified as inputs. It evaluates the matrix by panel
 * oriented assembly (\f$\ref{pc:ass}\f$) by first evaluating the interaction
 * matrix for all possible pairs of panels and then using the local to global
 * map of BEM spaces to fill the matrix entries.
 *
 * @param mesh ParametrizedMesh object containing all the parametrized
 *             panels in the mesh
 * @param space The trial and test BEM space to be used for evaluating
 *              the Galerkin matrix
 * @param N Order for Gauss Quadrature
 * @return An Eigen::MatrixXd type Galerkin Matrix for the given mesh and space
 */
Eigen::MatrixXd GalerkinMatrix(const ParametrizedMesh mesh,
                               const AbstractBEMSpace &space,
                               const unsigned int &N,
                               const double k);

/**
 * This function is used to evaluate the Single Layer Potential given by
 * \f$\Psi^{\Delta}_{SL}\Phi(x) = \int_{\Gamma} -\frac{1}{2\Pi} log ||x-y|| \Phi
 * (y) dS(y) \f$ for the function \f$\Phi = \sum_{i=1}^{N} c_{i} b^{i}_{N}\f$
 * where \f$c_{i}\f$ are the coefficients and \f$b^{i}_{N}\f$ are the basis
 * functions for the given BEM space and mesh. The Single Layer Potential is
 * evaluated using quadrature, at the evaluation point x passed as an input.
 *
 * @param x An Eigen::Vector2d type for the evaluation point
 * @param coeffs An Eigen::VectorXd type containing the coefficients \f$c_{i}\f$
 *               as mentioned above
 * @param mesh ParametrizedMesh object containing all the parametrized
 *             panels in the mesh
 * @param space The BEM space used for evaluating the Single Layer Potential
 * @param N Order for Gauss Quadrature
 * @return double representing the Single Layer Potential at the test point
 */
double Potential(const Eigen::Vector2d &x, const Eigen::VectorXd &coeffs,
                 const ParametrizedMesh &mesh, const AbstractBEMSpace &space,
                 const unsigned int &N);

} // namespace single_layer_correction
} // namespace parametricbem2d

#endif // SINGLELAYERHPP
