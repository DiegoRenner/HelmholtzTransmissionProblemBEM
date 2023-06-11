#include "mass_matrix.hpp"

namespace mass_matrix {

    Eigen::MatrixXd GalerkinMatrix(const ParametrizedMesh &mesh,
                                    const AbstractBEMSpace &trial_space,
                                    const AbstractBEMSpace &test_space,
                                    const QuadRule &GaussQR) {
        // Getting the number of panels in the mesh
        unsigned int numpanels = mesh.getNumPanels();
        // Getting the panels from the mesh
        const PanelVector &panels = mesh.getPanels();
        // Getting the number of local shape functions in the trial/test space
        unsigned int Qtrial = trial_space.getQ();
        unsigned int Qtest = test_space.getQ();
        // Getting dimensions of trial/test space
        unsigned int rows = test_space.getSpaceDim(numpanels);
        unsigned int cols = trial_space.getSpaceDim(numpanels);
        // Initializing the Galerkin matrix with zeros
        Eigen::MatrixXd output;
        output.setZero(rows, cols);
        // Panel oriented assembly
        std::vector<int> iv(std::min(rows, cols)), kv(GaussQR.n);
        std::iota(iv.begin(), iv.end(), 0);
        std::iota(kv.begin(), kv.end(), 0);
        std::for_each (iv.cbegin(), iv.cend(), [&](const int &i) {
            // Local to global mapping of the elements in interaction matrix
            auto &panel = *panels[i];
            size_t I, J, II, JJ;
            for (I = 0; I < Qtest; ++I) {
                for (J = 0; J < Qtrial; ++J) {
                    // Filling the Galerkin matrix entries
                    II = test_space.LocGlobMap(I + 1, i + 1, rows) - 1;
                    JJ = trial_space.LocGlobMap(J + 1, i + 1, cols) - 1;
                    output(II, JJ) += std::accumulate(kv.cbegin(), kv.cend(), 0., [&](const auto &sum, const int &k) {
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

} //namespace mass_matrix
