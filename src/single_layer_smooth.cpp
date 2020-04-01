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

#include "single_layer_smooth.hpp"
#include <math.h>

namespace parametricbem2d {
    namespace single_layer_smooth_helmholtz {

        typedef std::complex<double> complex_t;
        double sqrt_epsilon = sqrt(std::numeric_limits<double>::epsilon());
        double epsilon = std::numeric_limits<double>::epsilon();

        Eigen::MatrixXcd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                           const AbstractParametrizedCurve &pi_p,
                                           const AbstractBEMSpace &space,
                                           const QuadRule &GaussQR,
                                           const double k) {
            return ComputeIntegral(pi, pi_p, space, GaussQR, k);
        }

        Eigen::MatrixXcd ComputeIntegral(const AbstractParametrizedCurve &pi,
                                         const AbstractParametrizedCurve &pi_p,
                                         const AbstractBEMSpace &space,
                                         const QuadRule &GaussQR,
                                         const double k) {
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
                        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
                    };
                    auto G = [&](double s) { // Function associated with panel pi
                        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
                    };
                    auto integrand = [&](double s, double t) {
                        complex_t result = complex_t(0.0,0.0);
                        if ( (pi(s)-pi_p(t)).norm() > epsilon) {
                            result = complex_t(-y0(k*(pi(s)-pi_p(t)).norm())+(2/M_PI)*log(k*(pi(s)-pi_p(t)).norm()),
                                               j0(k*(pi(s)-pi_p(t)).norm()));
                        } else {
                            result = complex_t(0.0,
                                               j0(k*(pi(s)-pi_p(t)).norm()));
                        }
                        return result*F(t)*G(s);
                    };
                    complex_t integral = complex_t(0.,0.);
                    // Tensor product quadrature rule
                    for (unsigned int k = 0; k < N; ++k) {
                        for (unsigned int l = 0; l < N; ++l) {
                            double s = GaussQR.x(k);
                            double t = GaussQR.x(l);
                            double w = GaussQR.w(k)*GaussQR.w(l);
                            integral += w*integrand(s,t);
                        }
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = integral/4.;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &space,
                                        const unsigned int &N,
                                        const double k) {
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

    } //namespace single_layer_helmholtz
} // namespace parametricbem2d
