
#include <iostream>
#include "mass_matrix.hpp"

    namespace mass_matrix {

        typedef std::complex<double> complex_t;

        Eigen::MatrixXcd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                           const AbstractParametrizedCurve &pi_p,
                                           const AbstractBEMSpace &trial_space,
                                           const AbstractBEMSpace &test_space,
                                           const QuadRule &GaussQR,
                                           const QuadRule &CGaussQR) {
            if (&pi == &pi_p) { // Same Panels case
                return ComputeIntegral(pi, pi_p, trial_space, test_space, GaussQR);
            } else {
                unsigned int Qtrial = trial_space.getQ();
                unsigned int Qtest = test_space.getQ();
                // Interaction matrix with size Q x Q
                Eigen::MatrixXcd interaction_matrix = Eigen::MatrixXcd::Zero(Qtest,Qtrial);
                return interaction_matrix;
            }
        }


        Eigen::MatrixXcd ComputeIntegral(const AbstractParametrizedCurve &pi,
                                                const AbstractParametrizedCurve &pi_p,
                                                const AbstractBEMSpace &trial_space,
                                                const AbstractBEMSpace &test_space,
                                                const QuadRule &GaussQR){
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // Calculating the quadrature order for stable evaluation of integrands for
            // disjoint panels
            // No. of Reference Shape Functions in trial/test space
            int Qtrial = trial_space.getQ();
            int Qtest = test_space.getQ();
            // Interaction matrix with size Q x Q
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Computing the (i,j)th matrix entry
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    // Lambda expression for functions F and G
                    // Single Layer BIO
                    auto F = [&](double t) { // Function associated with panel pi_p
                        return trial_space.evaluateShapeFunction(j, t) * pi_p.Derivative_01(t).norm();
                    };
                    auto G = [&](double s) { // Function associated with panel pi
                        return test_space.evaluateShapeFunction(i, s);
                    };
                    complex_t integral = complex_t(0.,0.);
                    for (unsigned int k = 0; k < N; ++k) {
                        // Tensor product quadrature rule
                            double s = GaussQR.x(k);
                            integral += GaussQR.w(k)*F(s)*G(s);
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = integral;
                }
            }
            return interaction_matrix;
        }

        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                       const AbstractBEMSpace &trial_space,
                                        const AbstractBEMSpace &test_space,
                                       const unsigned int &N) {
                // Getting the number of panels in the mesh
            unsigned int numpanels = mesh.getNumPanels();
            // Getting dimensions of trial/test space
            // Getting the panels from the mesh
            PanelVector panels = mesh.getPanels();
            // Getting the number of local shape functions in the trial/test space
            unsigned int Qtrial = trial_space.getQ();
            unsigned int Qtest = test_space.getQ();
            unsigned int rows = test_space.getSpaceDim(numpanels);
            unsigned int cols = trial_space.getSpaceDim(numpanels);
            // Initializing the Galerkin matrix with zeros
            Eigen::MatrixXcd output = Eigen::MatrixXcd::Zero(rows, cols);
            // Panel oriented assembly
            QuadRule GaussQR = getGaussQR(N,0.,1.);
            QuadRule CGaussQR = getCGaussQR(N);
            for (unsigned int i = 0; i < numpanels; ++i) {
                for (unsigned int j = 0; j < numpanels; ++j) {
                    // Getting the interaction matrix for the pair of panels i and j
                    Eigen::MatrixXcd interaction_matrix =
                            InteractionMatrix(*panels[i], *panels[j], trial_space, test_space, GaussQR, CGaussQR);
                    // Local to global mapping of the elements in interaction matrix
                    for (unsigned int I = 0; I < Qtest; ++I) {
                        for (unsigned int J = 0; J < Qtrial; ++J) {
                            int II = test_space.LocGlobMap(I + 1, i + 1, numpanels) - 1;
                            int JJ = trial_space.LocGlobMap(J + 1, j + 1, numpanels) - 1;
                            // Filling the Galerkin matrix entries
                            output(II, JJ) += interaction_matrix(I, J);
                        }
                    }
                }
            }
            return output;
        }
    } //namespace single_layer_helmholtz
