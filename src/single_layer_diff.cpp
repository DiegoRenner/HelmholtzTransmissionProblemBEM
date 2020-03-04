/**
 * \file single_layer.cpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Single Layer BIO, using the transformations given in
 *        \f$\ref{ss:quadapprox}\f$ in the Lecture Notes for Advanced Numerical
 *        Methods for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "single_layer.hpp"

#include <limits>
#include <math.h>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"
#include "logweight_quadrature.hpp"
#include "parametrized_mesh.hpp"
#include <Eigen/Dense>
namespace parametricbem2d {
    namespace single_layer_helmholtz {

        typedef std::complex<double> complex_t;
        complex_t ii = complex_t(0,1.);
        double mascheroni = 0.57721566490153286060651; // More precision?
        double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        double tol = std::numeric_limits<double>::epsilon();
        //double sqrt_epsilon = std::numeric_limits<double>::epsilon();

        Eigen::MatrixXcd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                           const AbstractParametrizedCurve &pi_p,
                                           const AbstractBEMSpace &space,
                                           const QuadRule &GaussQR,
                                           const QuadRule &CGaussQR,
                                           const double c_o,
                                           const double c_i) {
            if (fabs(c_i) < tol){
                if (&pi == &pi_p) { // Same Panels case
                    return ComputeIntegralCoinciding(pi, pi_p, space, CGaussQR, c_o);
                }
                else if ((pi(1) - pi_p(-1)).norm() / 100. < tol ||
                         (pi(-1) - pi_p(1)).norm() / 100. < tol) {// Adjacent Panels case
                    return ComputeIntegralAdjacent(pi, pi_p, space, CGaussQR, c_o);
                }
                else {// Disjoint panels case
                    return ComputeIntegralGeneral(pi, pi_p, space, GaussQR, c_o);
                }
            } else {
                if (&pi == &pi_p) { // Same Panels case
                    return ComputeIntegralCoincidingDiff(pi, pi_p, space, CGaussQR, c_o, c_i);
                }
                    /*else if ((pi(1) - pi_p(-1)).norm() / 100. < tol ||
                             (pi(-1) - pi_p(1)).norm() / 100. < tol) {// Adjacent Panels case
                        return ComputeIntegralAdjacent(pi, pi_p, space, GaussQR, c_o);
                    }*/
                else {// Disjoint panels case
                    return ComputeIntegralGeneralDiff(pi, pi_p, space, GaussQR, c_o,c_i);
                }
            }
        }


        Eigen::MatrixXcd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                                 const AbstractParametrizedCurve &pi_p,
                                                 const AbstractBEMSpace &space,
                                                 const QuadRule &GaussQR,
                                                 const double c_o) {
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = space.getQ();
            // Interaction matrix with size Qtest x Qtrial
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Computing the (i,j)th matrix entry
            bool swap = (pi(1) - pi_p(-1)).norm() / 100. >
                        std::numeric_limits<double>::epsilon();
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:Kidp}\f$
                    auto F = [&](double t) {
                        if (swap) {
                            return space.evaluateShapeFunction_01_swapped(j, t) * pi_p.Derivative_01_swapped(t).norm();
                        } else {
                            return space.evaluateShapeFunction_01(j, t) * pi_p.Derivative_01(t).norm();
                        }
                    };

                    auto G = [&](double s) {
                        if (swap) {
                            return space.evaluateShapeFunction_01(i, s) * pi.Derivative_01(s).norm();
                        } else {
                            return space.evaluateShapeFunction_01_swapped(i, s) * pi.Derivative_01_swapped(s).norm();
                        }
                    };

                    // Lambda expression for the integrand in \f$\eqref{eq:Kidp}\f$
                    auto integrand = [&](double s, double t) {
                        complex_t k = complex_t(0.,0.);
                        // Finding the tangent of pi_p to get its normal

                        if (swap) {
                            if ( (pi[s]-pi_p.swapped_op(t)).norm() > tol) {
                                k = complex_t(-y0(c_o*(pi[s]-pi_p.swapped_op(t)).norm()),
                                              j0(c_o*(pi[s]-pi_p.swapped_op(t)).norm()));
                                //std::cout << "test" << std::endl;
                            };
                        } else {
                            if ( (pi.swapped_op(s)-pi_p[t]).norm() > tol) {
                                k = complex_t(-y0(c_o*(pi.swapped_op(s)-pi_p[t]).norm()),
                                              j0(c_o*(pi.swapped_op(s)-pi_p[t]).norm()));
                                //std::cout << "test" << std::endl;
                            };
                        }
                        return k * F(t) * G(s);
                    };

                    complex_t integral = complex_t(0.,0.);

                    // Tensor product quadrature for double integral in \f$\eqref{eq:Kidp}\f$
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            //double s = GaussQR.x(l)*(1.-GaussQR.x(k));
                            //double t = GaussQR.x(l);
                            //double w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
                            double s = GaussQR.x(k)*GaussQR.x(l);
                            double t = GaussQR.x(k);
                            double w = GaussQR.x(k)*GaussQR.w(k)*GaussQR.w(l);
                            integral += w*integrand(s,t);
                            double temp = t;
                            t = s;
                            s = temp;
                            integral += w*integrand(s,t);
                        }
                    }
                    // Filling the matrix entry
                    interaction_matrix(i, j) = 1/4.*integral;
                }
            }
            return interaction_matrix;
        }
        Eigen::MatrixXcd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                                   const AbstractParametrizedCurve &pi_p,
                                                   const AbstractBEMSpace &space,
                                                   const QuadRule &GaussQR,
                                                   const double c_o){
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // Calculating the quadrature order for stable evaluation of integrands for
            // disjoint panels as mentioned in \f$\ref{par:distpan}\f$
            // No. of Reference Shape Functions in trial/test space
            int Q = space.getQ();
            // Interaction matrix with size Q x Q
            Eigen::MatrixXcd interaction_matrix(Q, Q);
            // Computing the (i,j)th matrix entry

            for (int i = 0; i < Q; ++i) {
                for (int j = 0; j < Q; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
                    // Single Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunction_01(j, t) * pi_p.Derivative_01(t).norm();

                    };

                    auto G = [&](double s) { // Function associated with panel pi
                        return space.evaluateShapeFunction_01(i, s) * pi.Derivative_01(s).norm();
                    };

                    complex_t integral = complex_t(0.,0.);

                    // Tensor product quadrature rule
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            double s = GaussQR.x(l)*(1.-GaussQR.x(k));
                            double t = GaussQR.x(l);
                            double w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
                            if ( (pi[s]-pi_p[t]).norm() > tol) {
                                integral += complex_t(-y0(c_o*(pi[s]-pi_p[t]).norm()),
                                                      j0(c_o*(pi[s]-pi_p[t]).norm()))*F(t)*G(s)*w;
                            };
                            double temp = t;
                            t = s;
                            s = temp;
                            w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
                            if ( (pi[s]-pi_p[t]).norm() > tol) {
                                integral += complex_t(-y0(c_o*(pi[s]-pi_p[t]).norm()),
                                                      j0(c_o*(pi[s]-pi_p[t]).norm()))*F(t)*G(s)*w;
                            };
                        }
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = 1/4.*integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                                const AbstractParametrizedCurve &pi_p,
                                                const AbstractBEMSpace &space,
                                                const QuadRule &GaussQR,
                                                const double c_o){
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // Calculating the quadrature order for stable evaluation of integrands for
            // disjoint panels as mentioned in \f$\ref{par:distpan}\f$
            // No. of Reference Shape Functions in trial/test space
            int Q = space.getQ();
            // Interaction matrix with size Q x Q
            Eigen::MatrixXcd interaction_matrix(Q, Q);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Q; ++i) {
                for (int j = 0; j < Q; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
                    // Single Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunction_01(j, t) * pi_p.Derivative_01(t).norm();
                    };

                    auto G = [&](double s) { // Function associated with panel pi
                        //std::cout << pi.Derivative(s).norm() << std::endl;
                        return space.evaluateShapeFunction_01(i, s) * pi.Derivative_01(s).norm();
                    };

                    complex_t integral = complex_t(0.,0.);

                    // Tensor product quadrature rule
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            double s = GaussQR.x(k);
                            double t = GaussQR.x(l);
                            //integral += GaussQR.w(i) * GaussQR.w(j) *complex<double>(boost::math::sph_bessel(0,c_o*(pi(s)-pi_p(t)).norm()),
                            //        boost::math::sph_neumann(0,c_o*(pi(s)-pi_p(t)).norm())) * F(t) * G(s);
                            if ( (pi[s]-pi_p[t]).norm() > tol) {
                                integral += GaussQR.w(k)*GaussQR.w(l)*complex_t(-y0(c_o*(pi[s]-pi_p[t]).norm()),
                                                                                j0(c_o*(pi[s]-pi_p[t]).norm()))*F(t)*G(s);
                            };


                        }
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = 1/4.*integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd ComputeIntegralCoincidingDiff(const AbstractParametrizedCurve &pi,
                                                       const AbstractParametrizedCurve &pi_p,
                                                       const AbstractBEMSpace &space,
                                                       const QuadRule &GaussQR,
                                                       const double c_o,
                                                       const double c_i){
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // Calculating the quadrature order for stable evaluation of integrands for
            // disjoint panels as mentioned in \f$\ref{par:distpan}\f$
            // No. of Reference Shape Functions in trial/test space
            int Q = space.getQ();
            // Interaction matrix with size Q x Q
            Eigen::MatrixXcd interaction_matrix(Q, Q);
            // Computing the (i,j)th matrix entry

            for (int i = 0; i < Q; ++i) {
                for (int j = 0; j < Q; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
                    // Single Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunction_01(j, t) * pi_p.Derivative_01(t).norm();
                    };

                    auto G = [&](double s) { // Function associated with panel pi
                        //std::cout << pi.Derivative(s).norm() << std::endl;
                        return space.evaluateShapeFunction_01(i, s) * pi.Derivative_01(s).norm();
                    };

                    complex_t integral = complex_t(0.,0.);

                    // Tensor product quadrature rule
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            double s = GaussQR.x(l)*(1.-GaussQR.x(k));
                            double t = GaussQR.x(l);
                            double w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
                            if ((pi(s)-pi_p(t)).norm() > tol) {
                                integral += w * F(t) * G(s) * complex_t(-y0(c_o * (pi[s] - pi_p[t]).norm()) + y0(c_i * (pi[s] - pi_p[t]).norm()),
                                                                        j0(c_o * (pi[s] - pi_p[t]).norm()) - j0(c_i* (pi[s] - pi_p[t]).norm()));
                            }
                            double temp = t;
                            t = s;
                            s = temp;
                            w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
                            if ((pi(s)-pi_p(t)).norm() > tol) {
                                integral += w * F(t) * G(s) * complex_t(-y0(c_o * (pi[s] - pi_p[t]).norm()) + y0(c_i * (pi[s] - pi_p[t]).norm()),
                                                                        j0(c_o * (pi[s] - pi_p[t]).norm()) - j0(c_i* (pi[s] - pi_p[t]).norm()));
                            }
                        }
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = 1./4.*integral;
                }
            }
            return interaction_matrix;
        }
        Eigen::MatrixXcd ComputeIntegralGeneralDiff(const AbstractParametrizedCurve &pi,
                                                    const AbstractParametrizedCurve &pi_p,
                                                    const AbstractBEMSpace &space,
                                                    const QuadRule &GaussQR,
                                                    const double c_o,
                                                    const double c_i) {
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // No. of Reference Shape Functions in trial/test space
            int Q = space.getQ();
            // Interaction matrix with size Q x Q
            Eigen::MatrixXcd interaction_matrix(Q, Q);
            // Computing the (i,j)th matrix entry
            //std::cout << "test" << std::endl;
            for (int i = 0; i < Q; ++i) {
                for (int j = 0; j < Q; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
                    // Single Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
                    };

                    auto G = [&](double s) { // Function associated with panel pi
                        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
                    };

                    // Tensor product quadrature rule
                    complex_t integral = complex_t(0.,0.);
                    for (unsigned int i = 0; i < N; ++i) {
                        for (unsigned int j = 0; j < N; ++j) {
                            double s = GaussQR.x(i);
                            double t = GaussQR.x(j);
                            if ((pi(s)-pi_p(t)).norm() > tol) {
                                integral += GaussQR.w(i) * GaussQR.w(j) *
                                            (j0(c_o * (pi(s) - pi_p(t)).norm())
                                             + ii*y0(c_o * (pi(s) - pi_p(t)).norm()) -
                                             j0(c_i* (pi(s) - pi_p(t)).norm())
                                             - ii*y0(c_i * (pi(s) - pi_p(t)).norm())) *
                                            F(t) * G(s);
                            } else {
                                /*integral += GaussQR.w(i) * GaussQR.w(j) *
                                        (ii*1./(pi(s)-pi_p(t)).norm())*(1/c_o-1/c_i) *
                                           F(t) * G(s);*/
                            }
                        }
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = ii / 4. * integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &space,
                                        const unsigned int &N,
                                        const double c_o,
                                        const double c_i) {
            // Getting the number of panels in the mesh
            unsigned int numpanels = mesh.getNumPanels();
            // Getting dimensions of trial/test space
            unsigned int dims = space.getSpaceDim(numpanels);
            // Getting the panels from the mesh
            PanelVector panels = mesh.getPanels();
            // Getting the number of local shape functions in the trial/test space
            unsigned int Q = space.getQ();
            // Initializing the Galerkin matrix with zeros
            Eigen::MatrixXcd output = Eigen::MatrixXcd::Zero(dims, dims);
            // Panel oriented assembly \f$\ref{pc:ass}\f$
            QuadRule GaussQR = getGaussQR(N);
            QuadRule CGaussQR = getCGaussQR(N);
            for (unsigned int i = 0; i < numpanels; ++i) {
                for (unsigned int j = 0; j < numpanels; ++j) {
                    // std::cout << "For panels "<< i <<  "," <<j << " " ;
                    // Getting the interaction matrix for the pair of panels i and j
                    Eigen::MatrixXcd interaction_matrix =
                            InteractionMatrix(*panels[i], *panels[j], space, GaussQR, CGaussQR, c_o, c_i);
                    // Local to global mapping of the elements in interaction matrix
                    for (unsigned int I = 0; I < Q; ++I) {
                        for (unsigned int J = 0; J < Q; ++J) {
                            // int II = space.LocGlobMap(I + 1, i + 1, numpanels) - 1;
                            // int JJ = space.LocGlobMap(J + 1, j + 1, numpanels) - 1;
                            int II = space.LocGlobMap2(I + 1, i + 1, mesh) - 1;
                            int JJ = space.LocGlobMap2(J + 1, j + 1, mesh) - 1;
                            // Filling the Galerkin matrix entries
                            output(II, JJ) += interaction_matrix(I, J);
                        }
                    }
                }
            }
            return output;
        }

    } //namespace single_layer_helmholtz
} // namespace parametricbem2d
