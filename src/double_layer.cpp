/**
 * \file double_layer.cpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Double Layer BIO, using the transformations given in
 *        \f$\ref{ss:quadapprox}\f$ in the Lecture Notes for Advanced Numerical
 *        Methods for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package test
 */

#include "double_layer.hpp"

#include <limits>
#include <math.h>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "gauleg.hpp"
#include "logweight_quadrature.hpp"
#include "parametrized_mesh.hpp"
#include "/home/diegorenner/CLionProjects/Code/third_party/Betl2/Library/analytical_functions/hankel.hpp"

namespace parametricbem2d {
    namespace double_layer_helmholtz {
        typedef std::complex<double> complex_t;
        complex_t ii = complex_t(0,1);
        double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        double tol = std::numeric_limits<double>::epsilon();

        Eigen::MatrixXcd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                           const AbstractParametrizedCurve &pi_p,
                                           const AbstractBEMSpace &trial_space,
                                           const AbstractBEMSpace &test_space,
                                           const QuadRule &GaussQR,
                                           const QuadRule &CGaussQR,
                                           const double c_o,
                                           const double c_i) {

            if (fabs(c_i)<tol){
                if (&pi == &pi_p) { // Same Panels case
                    return ComputeIntegralCoinciding(pi, pi_p, trial_space, test_space, CGaussQR, c_o);
                }
                else if ((pi(1) - pi_p(-1)).norm() / 100. < tol ||
                         (pi(-1) - pi_p(1)).norm() / 100. < tol) {// Adjacent Panels case
                    return ComputeIntegralAdjacent(pi, pi_p, trial_space, test_space, CGaussQR, c_o);
                } else { //Disjoint panels case
                    return ComputeIntegralGeneral(pi, pi_p, trial_space, test_space, GaussQR, c_o);
                }
            } else {
                /*if (&pi == &pi_p) { // Same Panels case
                    return ComputeIntegralCoincidingDiff(pi, pi_p, trial_space, test_space, GaussQR, c_o, c_i);
                }
                else if ((pi(1) - pi_p(-1)).norm() / 100. < tol ||
                         (pi(-1) - pi_p(1)).norm() / 100. < tol) {// Adjacent Panels case
                    return ComputeIntegralAdjacentDiff(pi, pi_p, trial_space, test_space, GaussQR, c_o, c_i);
                }
                else {// Disjoint panels case*/
                return ComputeIntegralGeneralDiff(pi, pi_p, trial_space, test_space, CGaussQR, c_o, c_i);
                //}
            }
        }

        Eigen::MatrixXcd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                                   const AbstractParametrizedCurve &pi_p,
                                                   const AbstractBEMSpace &trial_space,
                                                   const AbstractBEMSpace &test_space,
                                                   const QuadRule &GaussQR,
                                                   const double c_o) {
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = test_space.getQ();
            // Interaction matrix with size Qtest x Qtrial
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:Kidp}\f$
                    auto F = [&](double t) {
                        return trial_space.evaluateShapeFunction_01(j, t) * pi_p.Derivative_01(t).norm();
                    };

                    auto G = [&](double s) {
                        return test_space.evaluateShapeFunction_01(i, s) * pi.Derivative_01(s).norm();
                    };

                    // Lambda expression for the integrand in \f$\eqref{eq:Kidp}\f$
                    auto integrand = [&](double s, double t) {
                        complex_t k = complex_t(0.,0.);
                        // Finding the tangent of pi_p to get its normal
                        Eigen::Vector2d tangent = pi_p.Derivative_01(t);
                        Eigen::Vector2d normal;
                        // Outward normal vector
                        normal << tangent(1), -tangent(0);
                        // Normalizing the normal vector
                        normal = normal / normal.norm();
                        if ( (pi[s]-pi_p[t]).norm() > tol && fabs((pi[s]-pi_p[t]).dot(normal)) > tol) { // Away from singularity
                            return (j1(c_o * (pi[s] - pi_p[t]).norm()) + ii * y1(c_o * (pi[s] - pi_p[t]).norm()))
                                   *c_o *((pi[s] - pi_p[t]).normalized()).dot(normal);
                        } else {
                            return complex_t(0.0,0.0);
                        }

                    };

                    complex_t integral = complex_t(0.,0.);

                    // Tensor product quadrature for double integral in \f$\eqref{eq:Kidp}\f$
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            double s = GaussQR.x(l)*(1.-GaussQR.x(k));
                            double t = GaussQR.x(l);
                            double w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
                            integral += w*integrand(s,t)*F(t)*G(s);
                            double temp = t;
                            t = s;
                            s = temp;
                            integral += w*integrand(s,t)*F(t)*G(s);
                        }
                    }
                    // Filling the matrix entry
                    interaction_matrix(i, j) = ii/4.*integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                                 const AbstractParametrizedCurve &pi_p,
                                                 const AbstractBEMSpace &trial_space,
                                                 const AbstractBEMSpace &test_space,
                                                 const QuadRule &GaussQR,
                                                 const double c_o) {
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = test_space.getQ();
            // Interaction matrix with size Qtest x Qtrial
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Computing the (i,j)th matrix entry
            bool swap = ((pi(1) - pi_p(-1)).norm() / 100. >
                        std::numeric_limits<double>::epsilon());
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:Kidp}\f$
                    auto F = [&](double t) {
                        if (swap) {
                            //double t_temp = swap ? 1-t:t;
                            return trial_space.evaluateShapeFunction_01_swapped(j, t) * pi_p.Derivative_01_swapped(t).norm();
                        } else {
                            return trial_space.evaluateShapeFunction_01(j, t) * pi_p.Derivative_01(t).norm();
                        }
                    };

                    auto G = [&](double s) {
                        //double s_temp = swap ? s:1-s;
                        if (swap) {
                            return test_space.evaluateShapeFunction_01(i, s) * pi.Derivative_01(s).norm();
                        } else {
                            return test_space.evaluateShapeFunction_01_swapped(i, s) * pi.Derivative_01_swapped(s).norm();
                        }
                    };

                    // Lambda expression for the integrand in \f$\eqref{eq:Kidp}\f$
                    auto integrand = [&](double s, double t) {
                        complex_t k = complex_t(0.,0.);
                        // Finding the tangent of pi_p to get its normal

                        Eigen::Vector2d tangent = swap ? pi_p.Derivative_01_swapped(t) : pi_p.Derivative_01(t);
                        Eigen::Vector2d normal;
                        // Outward normal vector
                        if (swap) {
                            normal << -tangent(1), tangent(0);
                        } else {
                            normal << tangent(1), -tangent(0);
                        }
                        // Normalizing the normal vector
                        normal = normal / normal.norm();
                        if (swap) {
                            if ( (pi[s]-pi_p.swapped_op(t)).norm() > tol && fabs((pi[s]-pi_p.swapped_op(t)).dot(normal)) > tol) { // Away from singularity
                                k = (j1(c_o * (pi[s] - pi_p.swapped_op(t)).norm()) +
                                     ii * y1(c_o * (pi[s] - pi_p.swapped_op(t)).norm()))
                                    * c_o *
                                    ((pi[s] - pi_p.swapped_op(t)).normalized()).dot(
                                            normal);
                            }
                        } else {
                            if ( (pi.swapped_op(s)-pi_p[t]).norm() > tol && fabs((pi.swapped_op(s)-pi_p[t]).dot(normal)) > tol) { // Away from singularity
                                k = (j1(c_o * (pi.swapped_op(s) - pi_p[t]).norm()) +
                                     ii * y1(c_o * (pi.swapped_op(s) - pi_p[t]).norm()))
                                    * c_o *
                                    ((pi.swapped_op(s) - pi_p[t]).normalized()).dot(
                                            normal);
                            }
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
                    interaction_matrix(i, j) = ii/4.*integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                                const AbstractParametrizedCurve &pi_p,
                                                const AbstractBEMSpace &trial_space,
                                                const AbstractBEMSpace &test_space,
                                                const QuadRule &GaussQR,
                                                const double c_o) {
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in space
            int Qtest = test_space.getQ();
            // Interaction matrix with size Qtest x Qtrial
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
                    // Double Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return trial_space.evaluateShapeFunction_01(j, t) *
                               pi_p.Derivative_01(t).norm();
                    };
                    auto G = [&](double s) { // Function associated with panel pi
                        return test_space.evaluateShapeFunction_01(i, s) *
                               pi.Derivative_01(s).norm();
                    };
                    // Lambda expression for \f$\hat{K}\f$ in \f$\eqref{eq:titg}\f$ for double
                    // Layer BIO
                    auto integrand = [&](double s, double t) {
                        // Finding the tangent of pi_p to get its normal
                        Eigen::Vector2d tangent = pi_p.Derivative_01(t);
                        Eigen::Vector2d normal;
                        // Outward normal vector
                        normal << tangent(1), -tangent(0);
                        // Normalizing the normal vector
                        normal = normal / normal.norm();
                        if ( (pi[s]-pi_p[t]).norm() > tol && fabs((pi[s]-pi_p[t]).dot(normal)) > tol) { // Away from singularity
                            return (j1(c_o * (pi[s] - pi_p[t]).norm()) + ii * y1(c_o * (pi[s] - pi_p[t]).norm()))
                                   *c_o *(pi[s] - pi_p[t]).normalized().dot(normal);

                        } else {
                            return complex_t(0., 0.);
                        }
                    };

                    complex_t integral = complex_t(0.,0.);

                    // Tensor product quadrature for double integral
                    for (unsigned int i = 0; i < N; ++i) {
                        for (unsigned int j = 0; j < N; ++j) {
                            double s = GaussQR.x(i);
                            double t = GaussQR.x(j);
                            integral +=
                                    GaussQR.w(i) * GaussQR.w(j) * integrand(s, t) * F(t) * G(s);
                        }
                    }
                    // Filling the matrix entry
                    interaction_matrix(i, j) = ii/4.*integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd ComputeIntegralGeneralDiff(const AbstractParametrizedCurve &pi,
                                                    const AbstractParametrizedCurve &pi_p,
                                                    const AbstractBEMSpace &trial_space,
                                                    const AbstractBEMSpace &test_space,
                                                    const QuadRule &GaussQR,
                                                    const double c_o,
                                                    const double c_i) {
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in space
            int Qtest = test_space.getQ();
            // Interaction matrix with size Qtest x Qtrial
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
                    // Double Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return trial_space.evaluateShapeFunction(j, t) *
                               pi_p.Derivative(t).norm();
                    };
                    auto G = [&](double s) { // Function associated with panel pi
                        return test_space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
                    };
                    // Lambda expression for \f$\hat{K}\f$ in \f$\eqref{eq:titg}\f$ for double
                    // Layer BIO
                    auto integrand = [&](double s, double t) {
                        // Finding the tangent of pi_p to get its normal
                        Eigen::Vector2d tangent = pi_p.Derivative(t);
                        Eigen::Vector2d normal;
                        // Outward normal vector
                        normal << tangent(1), -tangent(0);
                        // Normalizing the normal vector
                        normal = normal / normal.norm();
                        //if (fabs(c_o * (pi(s) - pi_p(t)).norm())>tol && fabs(c_i * (pi(s) - pi_p(t)).norm())>tol> tol) {
                        if ( (pi(s)-pi_p(t)).norm() > tol) {
                            return ((j1(c_o * (pi(s) - pi_p(t)).norm()) + ii*y1(c_o * (pi(s) - pi_p(t)).norm())) * c_o
                                    - (j1(c_i * (pi(s) - pi_p(t)).norm()) + ii*y1(c_i * (pi(s) - pi_p(t)).norm())) * c_i) *
                                   (pi(s) - pi_p(t)).normalized().dot(normal);
                        } else {
                            return complex_t(0.,0.);
                        };

                    };

                    complex_t integral = complex_t(0.,0.);

                    // Tensor product quadrature for double integral
                    for (unsigned int i = 0; i < N; ++i) {
                        for (unsigned int j = 0; j < N; ++j) {
                            double s = GaussQR.x(i);
                            double t = GaussQR.x(j);
                            integral +=
                                    GaussQR.w(i) * GaussQR.w(j) * integrand(s, t) * F(t) * G(s);
                        }
                    }
                    // Filling the matrix entry
                    interaction_matrix(i, j) = ii/4.*integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &trial_space,
                                        const AbstractBEMSpace &test_space,
                                        const unsigned int &N,
                                        const double c_o,
                                        const double c_i) {
            // Getting number of panels in the mesh
            unsigned int numpanels = mesh.getNumPanels();
            // Getting dimensions for trial and test spaces
            unsigned int rows = test_space.getSpaceDim(numpanels);
            unsigned int cols = trial_space.getSpaceDim(numpanels);
            // Getting the panels from the mesh
            PanelVector panels = mesh.getPanels();
            // Getting the number of local shape functions in the trial and test spaces
            unsigned int Qtest = test_space.getQ();
            unsigned int Qtrial = trial_space.getQ();
            // Initializing the Galerkin matrix with zeros
            Eigen::MatrixXcd output = Eigen::MatrixXd::Zero(rows, cols);
            // Panel oriented assembly \f$\ref{pc:ass}\f$
            //QuadRule LogWeightQR = getLogWeightQR(1, N);
            QuadRule GaussQR = getGaussQR(N);
            QuadRule CGaussQR = getCGaussQR(N);
            for (unsigned int i = 0; i < numpanels; ++i) {
                for (unsigned int j = 0; j < numpanels; ++j) {
                    // Getting the interaction matrix for the pair of panels i and j
                    Eigen::MatrixXcd interaction_matrix = InteractionMatrix(
                            *panels[i], *panels[j], trial_space, test_space, GaussQR, CGaussQR, c_o, c_i);
                    // Local to global mapping of the elements in interaction matrix
                    for (unsigned int I = 0; I < Qtest; ++I) {
                        for (unsigned int J = 0; J < Qtrial; ++J) {
                            //int II = test_space.LocGlobMap(I + 1, i + 1, numpanels) - 1;
                            //int JJ = trial_space.LocGlobMap(J + 1, j + 1, numpanels) - 1;
                            int II = test_space.LocGlobMap2(I + 1, i + 1, mesh) - 1;
                            int JJ = trial_space.LocGlobMap2(J + 1, j + 1, mesh) - 1;
                            // Filling the Galerkin matrix entries
                            output(II, JJ) += interaction_matrix(I, J);
                        }
                    }
                }
            }
            return output;
        }

    }
} // namespace parametricbem2d

