#include "mass_matrix.hpp"
#include <execution>

    namespace mass_matrix {

        typedef std::complex<double> complex_t;
        static const complex_t czero(0., 0.);
#if 0
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
                return Eigen::MatrixXcd::Zero(Qtest, Qtrial);
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
                    complex_t integral = czero;
                    for (unsigned int k = 0; k < N; ++k) {
                        // Tensor product quadrature rule
                        double s = GaussQR.x(k);
                        auto F = trial_space.evaluateShapeFunction(j, s) * pi_p.Derivative_01(s).norm();
                        auto G = test_space.evaluateShapeFunction(i, s);
                        integral += GaussQR.w(k) * F * G;
                    }
                    // Filling up the matrix entry
                    interaction_matrix(i, j) = integral;
                }
            }
            return interaction_matrix;
        }
#endif
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
            QuadRule GaussQR = getGaussQR(N, 0., 1.);
            QuadRule CGaussQR = getCGaussQR(N);
            std::vector<int> iv(numpanels), kv(GaussQR.n);
            std::iota(iv.begin(), iv.end(), 0);
            std::iota(kv.begin(), kv.end(), 0);
            std::for_each (std::execution::par_unseq, iv.cbegin(), iv.cend(), [&](const int &i) {
                // Local to global mapping of the elements in interaction matrix
                auto &panel = *panels[i];
                for (unsigned int I = 0; I < Qtest; ++I) {
                    for (unsigned int J = 0; J < Qtrial; ++J) {
                        // Filling the Galerkin matrix entries
                        output(test_space.LocGlobMap(I + 1, i + 1, numpanels) - 1, trial_space.LocGlobMap(J + 1, i + 1, numpanels) - 1) +=
                        std::accumulate(kv.cbegin(), kv.cend(), czero, [&](const auto &sum, const int &k) {
                            // Tensor product quadrature rule
                            double s = GaussQR.x(k);
                            auto F = trial_space.evaluateShapeFunction(J, s) * panel.Derivative_01(s).norm();
                            auto G = test_space.evaluateShapeFunction(I, s);
                            return sum + GaussQR.w(k) * F * G;
                        });
                    }
                }
            });
            return output;
        }
    } //namespace single_layer_helmholtz
