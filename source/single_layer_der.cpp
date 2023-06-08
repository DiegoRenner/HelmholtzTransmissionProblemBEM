#include "single_layer_der.hpp"
#include "cbessel.hpp"
#include <execution>

namespace single_layer_helmholtz_der {

    typedef std::complex<double> complex_t;
    static const complex_t ii(0., 1.);
    static const complex_t czero(0., 0.);
    static const double epsilon = std::numeric_limits<double>::epsilon();

    inline
    void ComputeIntegralAdjacent(Eigen::MatrixXcd &interaction_matrix,
                                 const AbstractParametrizedCurve &pi,
                                 const AbstractParametrizedCurve &pi_p,
                                 const AbstractBEMSpace &space,
                                 const QuadRule &GaussQR,
                                 const complex_t &k,
                                 double c_i, double c_o,
                                 gq_workspace_t &ws) throw() {
        unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
        // The number of Reference Shape Functions in trial space
        int Qtrial = space.getQ();
        // The number of Reference Shape Functions in test space
        int Qtest = space.getQ();
        bool swap = (pi(1) - pi_p(-1)).norm() / 100. > epsilon;
        double sqrtc_i = sqrt(c_i), sqrtc_o = sqrt(c_o);
        complex_t ksqrtc_i = k * sqrtc_i, ksqrtc_o = k * sqrtc_o;
        double ksqrtca_i = std::abs(ksqrtc_i), ksqrtca_o = std::abs(ksqrtc_o);
        unsigned N2 = N * N;
        bool with_i = c_i != 0.;
        // compute the required values once
        std::for_each (std::execution::par_unseq, ws.ind().cbegin(), ws.ind().cbegin() + N2, [&](auto &i) {
            unsigned I = i / N, J = i - I * N;
            ws.IJ(i) = std::make_pair(I, J);
            complex_t &result_i = ws.result_i(i), &result_o = ws.result_o(i);
            result_i = result_o = czero;
            double &t = ws.t(i), &s = ws.s(i), &w = ws.w(i), &t_norm = ws.t_norm(i), &tp_norm = ws.tp_norm(i);
            t = GaussQR.x(I);
            s = t * GaussQR.x(J);
            w = t * GaussQR.w(I) * GaussQR.w(J);
            tp_norm = swap ? pi_p.Derivative_01_swapped(t).norm() : pi_p.Derivative_01(t).norm();
            t_norm = swap ? pi.Derivative_01(s).norm() : pi.Derivative_01_swapped(s).norm();
            double d = swap ? (pi[s]-pi_p.swapped_op(t)).norm() : (pi.swapped_op(s)-pi_p[t]).norm();
            if (with_i && ksqrtca_i * d > epsilon)
                result_i = complex_bessel::H1(1, ksqrtc_i * d) * d;
            if (ksqrtca_o * d > epsilon)
                result_o = complex_bessel::H1(1, ksqrtc_o * d) * d;
        });
        // Lambda expression for the integrand
        auto integrand = [&](int i, int j, int m) {
            const double &s = ws.s(m), &t = ws.t(m);
            auto F = (swap ? space.evaluateShapeFunction_01_swapped(j, t) : space.evaluateShapeFunction(j, t)) * ws.tp_norm(m);
            auto G = (swap ? space.evaluateShapeFunction(i, s) : space.evaluateShapeFunction_01_swapped(i, s)) * ws.t_norm(m);
            return (sqrtc_o * ws.result_o(m) - sqrtc_i * ws.result_i(m)) * F * G;
        };
        for (int i = 0; i < Qtest; ++i) {
            for (int j = 0; j < Qtrial; ++j) {
                interaction_matrix(i, j) =
                -0.25 * ii * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                    const int &I = ij.first, &J = ij.second;
                    return sum + ws.w(I * N + J) * (integrand(i, j, I * N + J) + integrand(i, j, J * N + I));
                });
            }
        }
    }

    inline
    void ComputeIntegralCoinciding(Eigen::MatrixXcd &interaction_matrix,
                                   const AbstractParametrizedCurve &pi,
                                   const AbstractParametrizedCurve &pi_p,
                                   const AbstractBEMSpace &space,
                                   const QuadRule &GaussQR,
                                   const complex_t &k,
                                   double c_i, double c_o,
                                   gq_workspace_t &ws) throw() {
        unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
        // Calculating the quadrature order for stable evaluation of integrands for
        // disjoint panels
        // No. of Reference Shape Functions in trial/test space
        int Q = space.getQ();
        double sqrtc_i = sqrt(c_i), sqrtc_o = sqrt(c_o);
        complex_t ksqrtc_i = k * sqrtc_i, ksqrtc_o = k * sqrtc_o;
        double ksqrtca_i = std::abs(ksqrtc_i), ksqrtca_o = std::abs(ksqrtc_o);
        unsigned N2 = N * N;
        bool with_i = c_i != 0.;
        // compute the required values once
        std::for_each (std::execution::par_unseq, ws.ind().cbegin(), ws.ind().cbegin() + N2, [&](auto &i) {
            unsigned I = i / N, J = i - I * N;
            ws.IJ(i) = std::make_pair(I, J);
            complex_t &result_i = ws.result_i(i), &result_o = ws.result_o(i);
            result_i = result_o = czero;
            double &t = ws.t(i), &s = ws.s(i), &w = ws.w(i), &t_norm = ws.t_norm(i), &tp_norm = ws.tp_norm(i);
            t = GaussQR.x(J);
            s = t * (1. - GaussQR.x(I));
            w = t * GaussQR.w(I) * GaussQR.w(J);
            tp_norm = pi_p.Derivative_01(t).norm();
            t_norm = pi.Derivative_01(s).norm();
            double d = (pi[s] - pi_p[t]).norm();
            if (with_i && ksqrtca_i * d > epsilon)
                result_i = complex_bessel::H1(1, ksqrtc_i * d) * d;
            if (ksqrtca_o * d > epsilon)
                result_o = complex_bessel::H1(1, ksqrtc_o * d) * d;
        });
        auto integrand = [&](int i, int j, int m) {
            auto F = space.evaluateShapeFunction(j, ws.t(m)) * ws.tp_norm(m);
            auto G = space.evaluateShapeFunction(i, ws.s(m)) * ws.t_norm(m);
            return (sqrtc_o * ws.result_o(m) - sqrtc_i * ws.result_i(m)) * F * G;
        };
        for (int i = 0; i < Q; ++i) {
            for (int j = 0; j < Q; ++j) {
                interaction_matrix(i, j) =
                -0.25 * ii * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                    const int &I = ij.first, &J = ij.second;
                    return sum + ws.w(I * N + J) * (integrand(i, j, I * N + J) + integrand(i, j, J * N + I));
                });
            }
        }
    }

    inline
    void ComputeIntegralGeneral(Eigen::MatrixXcd &interaction_matrix,
                                const AbstractParametrizedCurve &pi,
                                const AbstractParametrizedCurve &pi_p,
                                const AbstractBEMSpace &space,
                                const QuadRule &GaussQR,
                                const complex_t &k,
                                double c_i, double c_o,
                                gq_workspace_t &ws) throw() {
        unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
        // Calculating the quadrature order for stable evaluation of integrands for
        // disjoint panels
        // No. of Reference Shape Functions in trial/test space
        int Q = space.getQ();
        double sqrtc_i = sqrt(c_i), sqrtc_o = sqrt(c_o);
        complex_t ksqrtc_i = k * sqrtc_i, ksqrtc_o = k * sqrtc_o;
        double ksqrtca_i = std::abs(ksqrtc_i), ksqrtca_o = std::abs(ksqrtc_o);
        unsigned N2 = N * N;
        bool with_i = c_i != 0.;
        // compute the required values once
        std::for_each (std::execution::par_unseq, ws.ind().cbegin(), ws.ind().cbegin() + N2, [&](auto &i) {
            unsigned I = i / N, J = i - I * N;
            ws.IJ(i) = std::make_pair(I, J);
            complex_t &result_i = ws.result_i(i), &result_o = ws.result_o(i);
            result_i = result_o = czero;
            double &t = ws.t(i), &s = ws.s(i), &w = ws.w(i), &t_norm = ws.t_norm(i), &tp_norm = ws.tp_norm(i);
            s = GaussQR.x(I);
            t = GaussQR.x(J);
            w = GaussQR.w(I) * GaussQR.w(J);
            tp_norm = pi_p.Derivative_01(t).norm();
            t_norm = pi.Derivative_01(s).norm();
            double d = (pi[s] - pi_p[t]).norm();
            if (with_i && ksqrtca_i * d > epsilon)
                result_i = complex_bessel::H1(1, ksqrtc_i * d) * d;
            if (ksqrtca_o * d > epsilon)
                result_o = complex_bessel::H1(1, ksqrtc_o * d) * d;
        });
        for (int i = 0; i < Q; ++i) {
            for (int j = 0; j < Q; ++j) {
                // Filling up the matrix entry
                interaction_matrix(i, j) =
                -0.25 * ii * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                    int m = ij.first * N + ij.second;
                    auto F = space.evaluateShapeFunction(j, ws.t(m)) * ws.tp_norm(m);
                    auto G = space.evaluateShapeFunction(i, ws.s(m)) * ws.t_norm(m);
                    return sum + ws.w(m) * (sqrtc_o * ws.result_o(m) - sqrtc_i * ws.result_i(m)) * F * G;
                });
            }
        }
    }

    inline
    void InteractionMatrix(Eigen::MatrixXcd &interaction_matrix,
                           const AbstractParametrizedCurve &pi,
                           const AbstractParametrizedCurve &pi_p,
                           const AbstractBEMSpace &space,
                           const QuadRule &GaussQR,
                           const QuadRule &CGaussQR,
                           const complex_t &k,
                           double c_i, double c_o,
                           gq_workspace_t &ws) throw() {
        if (&pi == &pi_p) { // Same Panels case
            ComputeIntegralCoinciding(interaction_matrix, pi, pi_p, space, CGaussQR, k, c_i, c_o, ws);
        }
        else if ((pi(1) - pi_p(-1)).norm() / 100. < epsilon ||
                    (pi(-1) - pi_p(1)).norm() / 100. < epsilon) {// Adjacent Panels case
            ComputeIntegralAdjacent(interaction_matrix, pi, pi_p, space, CGaussQR, k, c_i, c_o, ws);
        }
        else {// Disjoint panels case
            ComputeIntegralGeneral(interaction_matrix, pi, pi_p, space, GaussQR, k, c_i, c_o, ws);
        }
    }

    Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh &mesh,
                                    const AbstractBEMSpace &space,
                                    const QuadRule &GaussQR,
                                    const QuadRule &CGaussQR,
                                    const complex_t &k,
                                    double c_i,
                                    double c_o) {
        // Getting the number of panels in the mesh
        unsigned int numpanels = mesh.getNumPanels();
        // Getting dimensions of trial/test space
        unsigned int dims = space.getSpaceDim(numpanels);
        // Getting the panels from the mesh
        const PanelVector &panels = mesh.getPanels();
        // Getting the number of local shape functions in the trial/test space
        unsigned int Q = space.getQ();
        Eigen::MatrixXcd interaction_matrix(Q, Q);
        // Initializing the Galerkin matrix with zeros
        Eigen::MatrixXcd output = Eigen::MatrixXcd::Zero(dims, dims);
        // Panel oriented assembly
        gq_workspace_t ws(std::max(GaussQR.n, CGaussQR.n));
        unsigned i, j, I, J;
        for (i = 0; i < numpanels; ++i) {
            const auto &pi = *panels[i];
            for (j = 0; j < numpanels; ++j) {
                const auto &pj = *panels[j];
                // Getting the interaction matrix for the pair of panels i and j
                InteractionMatrix(interaction_matrix, pi, pj, space, GaussQR, CGaussQR, k, c_i, c_o, ws);
                // Local to global mapping of the elements in interaction matrix
                for (I = 0; I < Q; ++I) {
                    for (J = 0; J < Q; ++J) {
                        // Filling the Galerkin matrix entries
                        output(space.LocGlobMap(I + 1, i + 1, numpanels) - 1, space.LocGlobMap(J + 1, j + 1, numpanels) - 1)
                            += interaction_matrix(I, J);
                    }
                }
            }
        }
        return output;
    }

} //namespace single_layer_helmholtz
