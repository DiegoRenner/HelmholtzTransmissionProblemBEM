/**
 * \file hypersingular.cpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Hypersingular BIO, using the transformations given in
 *        \f$\ref{ss:quadapprox}\f$ in the Lecture Notes for Advanced Numerical
 *        Methods for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "hypersingular_smooth.hpp"

namespace parametricbem2d {
    namespace hypersingular_smooth_helmholtz {

        typedef std::complex<double> complex_t;
        double epsilon = std::numeric_limits<double>::epsilon();

        Eigen::MatrixXcd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                           const AbstractParametrizedCurve &pi_p,
                                           const AbstractBEMSpace &space,
                                           const QuadRule &GaussQR,
                                           const double k){
            return ComputeIntegralGeneral(pi, pi_p, space, GaussQR, k);
        }

        Eigen::MatrixXcd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                                const AbstractParametrizedCurve &pi_p,
                                                const AbstractBEMSpace &space,
                                                const QuadRule &GaussQR,
                                                const double k) {
            unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in space
            int Q = space.getQ();
            // The number of Reference Shape Functions in space
            // Interaction matrix with size Qtest x Qtrial
            Eigen::MatrixXcd interaction_matrix(Q, Q);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Q; ++i) {
                for (int j = 0; j < Q; ++j) {
                    // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
                    // Double Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunction_01(j, t)*pi_p.Derivative_01(t).norm();
                    };
                    auto G = [&](double s) { // Function associated with panel pi
                        return space.evaluateShapeFunction_01(i, s)*pi.Derivative_01(s).norm();
                    };
                    auto F_arc = [&](double t) { // Function associated with panel pi_p
                        return space.evaluateShapeFunctionDot_01(j, t);
                    };
                    auto G_arc = [&](double s) { // Function associated with panel pi
                        return space.evaluateShapeFunctionDot_01(i, s);
                    };
                    // Lambda expression for \f$\hat{K}\f$ in \f$\eqref{eq:titg}\f$ for double
                    // Layer BIO
                    auto integrand = [&](double s, double t) {
                        complex_t result = complex_t(0.,0.);
                        // Finding the tangent of pi_p to get its normal
                        Eigen::Vector2d tangent_p = pi_p.Derivative_01(t);
                        Eigen::Vector2d normal_p;
                        // Outward normal vector
                        normal_p << tangent_p(1), -tangent_p(0);
                        // Normalizing the normal vector
                        normal_p = normal_p / normal_p.norm();
                        Eigen::Vector2d tangent = pi.Derivative_01(s);
                        Eigen::Vector2d normal;
                        // Outward normal vector
                        normal << tangent(1), -tangent(0);
                        // Normalizing the normal vector
                        normal = normal / normal.norm();
                        if ( (pi(s)-pi_p(t)).norm() > epsilon) {
                            result = complex_t(-y0(k*(pi(s)-pi_p(t)).norm())+(2/M_PI)*log(k*(pi(s)-pi_p(t)).norm()),
                                               j0(k*(pi(s)-pi_p(t)).norm()));
                        } else {
                            result = complex_t(0.0,
                                               j0(k*(pi(s)-pi_p(t)).norm()));
                        }
                        return result * (F_arc(t) * G_arc(s) - k * k * F(t) * G(s) * normal.dot(normal_p));
                    };
                    complex_t integral = complex_t(0.,0.);
                    // Tensor product quadrature for double integral
                    for (unsigned int i = 0; i < N; ++i) {
                        for (unsigned int j = 0; j < N; ++j) {
                            double s = GaussQR.x(i);
                            double t = GaussQR.x(j);
                            double w = GaussQR.w(i)*GaussQR.w(j);
                            integral += w*integrand(s, t);
                        }
                    }
                    // Filling the matrix entry
                    interaction_matrix(i, j) = integral/4.;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &space,
                                        const unsigned int &N,
                                        const double k){
            // Getting the number of panels in the mesh
            unsigned int numpanels = mesh.getNumPanels();
            // Getting dimensions of trial/test space
            unsigned int dims = space.getSpaceDim(numpanels);
            // Getting the panels from the mesh
            PanelVector panels = mesh.getPanels();
            // Getting the number of local shape functions in the trial/test space
            unsigned int Q = space.getQ();
            // Initializing the Galerkin matrix with zeros
            Eigen::MatrixXcd output = Eigen::MatrixXd::Zero(dims, dims);
            // Panel oriented assembly \f$\ref{pc:ass}\f$
            //QuadRule LogWeightQR = getLogWeightQR(1, N);
            QuadRule GaussQR = getGaussQR(N,-1.,1.);
            for (unsigned int i = 0; i < numpanels; ++i) {
                for (unsigned int j = 0; j < numpanels; ++j) {
                    // Getting the interaction matrix for the pair of panels i and j
                    Eigen::MatrixXcd interaction_matrix =
                            InteractionMatrix(*panels[i], *panels[j], space, GaussQR, k);
                    // Local to global mapping of the elements in interaction matrix
                    for (unsigned int I = 0; I < Q; ++I) {
                        for (unsigned int J = 0; J < Q; ++J) {
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

    }// namespace hypersingular_helmholtz
} // namespace parametricbem2d
