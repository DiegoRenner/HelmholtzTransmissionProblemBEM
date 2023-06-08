#include "mass_matrix.hpp"
#include <execution>

namespace mass_matrix {

    typedef std::complex<double> complex_t;
    static const complex_t czero(0., 0.);

    Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh &mesh,
                                    const AbstractBEMSpace &trial_space,
                                    const AbstractBEMSpace &test_space,
                                    const QuadRule &GaussQR) {
        // Getting the number of panels in the mesh
        unsigned int numpanels = mesh.getNumPanels();
        // Getting dimensions of trial/test space
        // Getting the panels from the mesh
        const PanelVector &panels = mesh.getPanels();
        // Getting the number of local shape functions in the trial/test space
        unsigned int Qtrial = trial_space.getQ();
        unsigned int Qtest = test_space.getQ();
        unsigned int rows = test_space.getSpaceDim(numpanels);
        unsigned int cols = trial_space.getSpaceDim(numpanels);
        // Initializing the Galerkin matrix with zeros
        Eigen::MatrixXcd output = Eigen::MatrixXcd::Zero(rows, cols);
        // Panel oriented assembly
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
