#include "hypersingular_der2.hpp"
#include "discontinuous_space.hpp"
#include "cbessel.hpp"
#include <execution>

namespace hypersingular_helmholtz_der2 {

    typedef std::complex<double> complex_t;
    static const complex_t ii(0., 1.);
    static const complex_t czero(0., 0.);
    static const double epsilon = std::numeric_limits<double>::epsilon();

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
        // The number of Reference Shape Functions in trial space
        int Q = space.getQ();
        complex_t ksqrtc_i = k * sqrt(c_i), ksqrtc_o = k * sqrt(c_o), k2c_i = k * k * c_i, k2c_o = k * k * c_o;
        double ksqrtca_i = std::abs(ksqrtc_i), ksqrtca_o = std::abs(ksqrtc_o);
        unsigned N2 = N * N;
        bool with_i = c_i != 0.;
        // compute the required values once
        std::for_each (std::execution::par_unseq, ws.ind().cbegin(), ws.ind().cbegin() + N2, [&](auto &i) {
            unsigned I = i / N, J = i - I * N;
            ws.IJ(i) = std::make_pair(I, J);
            complex_t &result_i = ws.result_i(i), &result_o = ws.result_o(i), &result_ir = ws.result_i(i, true), &result_or = ws.result_o(i, true);
            result_i = result_o = result_ir = result_or = czero;
            double &t = ws.t(i), &s = ws.s(i), &w = ws.w(i), &t_norm = ws.t_norm(i), &tp_norm = ws.tp_norm(i), &n_dot_np = ws.n_dot_np(i);
            t = GaussQR.x(J);
            s = t * (1. - GaussQR.x(I));
            w = t * GaussQR.w(I) * GaussQR.w(J);
            // Finding the tangent of pi_p to get its normal
            const Eigen::Vector2d &tangent_p = pi_p.Derivative_01(t);
            const Eigen::Vector2d &tangent = pi.Derivative_01(s);
            tp_norm = tangent_p.norm();
            t_norm = tangent.norm();
            n_dot_np = (tangent_p(1) * tangent(1) + tangent_p(0) * tangent(0)) / (tp_norm * t_norm);
            double d = (pi[s] - pi_p[t]).norm();
            if (with_i && ksqrtca_i * d > epsilon) {
                complex_t h0 = complex_bessel::H1(0, ksqrtc_i * d), h1 = complex_bessel::H1(1, ksqrtc_i * d);
                result_i = h0 * d * d - h1 * d / ksqrtc_i;
                result_ir = -h1 * d * 4.0 * ksqrtc_i + h0 * 2.0;
            }
            if (ksqrtca_o * d > epsilon) {
                complex_t h0 = complex_bessel::H1(0, ksqrtc_o * d), h1 = complex_bessel::H1(1, ksqrtc_o * d);
                result_o = h0 * d * d - h1 * d / ksqrtc_o;
                result_or = -h1 * d * 4.0 * ksqrtc_o + h0 * 2.0;
            }
        });
        // Lambda expression for the integrand
        auto integrand = [&](int i, int j, int m) {
            const double &s = ws.s(m), &t = ws.t(m);
            auto F = space.evaluateShapeFunction(j, t) * ws.tp_norm(m);
            auto G = space.evaluateShapeFunction(i, s) * ws.t_norm(m);
            auto F_arc = space.evaluateShapeFunctionDot_01(j, t);
            auto G_arc = space.evaluateShapeFunctionDot_01(i, s);
            complex_t r_o1 = c_o * ws.result_o(m), r_o2 = c_o * ws.result_o(m, true);
            complex_t r_i1 = c_i * ws.result_i(m), r_i2 = c_i * ws.result_i(m, true);
            return (r_o1 - r_i1) * F_arc * G_arc + (r_o2 - r_i2 - k2c_o * r_o1 + k2c_i * r_i1) * F * G * ws.n_dot_np(m);
        };
        for (int i = 0; i < Q; ++i) {
            for (int j = 0; j < Q; ++j) {
                // Filling the matrix entry
                interaction_matrix(i, j) =
                -ii * 0.25 * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                    const int &I = ij.first, &J = ij.second;
                    return sum + ws.w(I * N + J) * (integrand(i, j, I * N + J) + integrand(i, j, J * N + I));
                });
            }
        }
    }

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
        bool swap = ((pi(1) - pi_p(-1)).norm() / 100. > epsilon);
        complex_t ksqrtc_i = k * sqrt(c_i), ksqrtc_o = k * sqrt(c_o), k2c_i = k * k * c_i, k2c_o = k * k * c_o;
        double ksqrtca_i = std::abs(ksqrtc_i), ksqrtca_o = std::abs(ksqrtc_o);
        unsigned N2 = N * N;
        bool with_i = c_i != 0.;
        // compute the required values once
        std::for_each (std::execution::par_unseq, ws.ind().cbegin(), ws.ind().cbegin() + N2, [&](auto &i) {
            unsigned I = i / N, J = i - I * N;
            ws.IJ(i) = std::make_pair(I, J);
            complex_t &result_i = ws.result_i(i), &result_o = ws.result_o(i), &result_ir = ws.result_i(i, true), &result_or = ws.result_o(i, true);
            result_i = result_o = result_ir = result_or = czero;
            double &t = ws.t(i), &s = ws.s(i), &w = ws.w(i), &t_norm = ws.t_norm(i), &tp_norm = ws.tp_norm(i), &n_dot_np = ws.n_dot_np(i);
            t = GaussQR.x(I);
            s = t * GaussQR.x(J);
            w = t * GaussQR.w(I) * GaussQR.w(J);
            // Finding the tangent of pi_p to get its normal
            const Eigen::Vector2d &tangent_p = swap ? -pi_p.Derivative_01_swapped(t) : pi_p.Derivative_01(t);
            const Eigen::Vector2d &tangent = swap ? pi.Derivative_01(s) : -pi.Derivative_01_swapped(s);
            tp_norm = tangent_p.norm();
            t_norm = tangent.norm();
            n_dot_np = (tangent_p(1) * tangent(1) + tangent_p(0) * tangent(0)) / (tp_norm * t_norm);
            double d = swap ? (pi[s] - pi_p.swapped_op(t)).norm() : (pi.swapped_op(s) - pi_p[t]).norm();
            if (with_i && ksqrtca_i * d > epsilon) {
                complex_t h0 = complex_bessel::H1(0, ksqrtc_i * d), h1 = complex_bessel::H1(1, ksqrtc_i * d);
                result_i = h0 * d * d - h1 * d / ksqrtc_i;
                result_ir = -h1 * d * 4.0 * ksqrtc_i + h0 * 2.0;
            }
            if (ksqrtca_o * d > epsilon) {
                complex_t h0 = complex_bessel::H1(0, ksqrtc_o * d), h1 = complex_bessel::H1(1, ksqrtc_o * d);
                result_o = h0 * d * d - h1 * d / ksqrtc_o;
                result_or = -h1 * d * 4.0 * ksqrtc_o + h0 * 2.0;
            }
        });
        // Lambda expression for the integrand
        auto integrand = [&](int i, int j, int m) {
            const double &s = ws.s(m), &t = ws.t(m);
            auto F = (swap ? space.evaluateShapeFunction_01_swapped(j, t) : space.evaluateShapeFunction(j, t)) * ws.tp_norm(m);
            auto G = (swap ? space.evaluateShapeFunction(i, s) : space.evaluateShapeFunction_01_swapped(i, s)) * ws.t_norm(m);
            auto F_arc = swap ? space.evaluateShapeFunctionDot_01_swapped(j, t) : space.evaluateShapeFunctionDot_01(j, t);
            auto G_arc = swap ? space.evaluateShapeFunctionDot_01(i, s) : space.evaluateShapeFunctionDot_01_swapped(i, s);
            complex_t r_o1 = c_o * ws.result_o(m), r_o2 = c_o * ws.result_o(m, true);
            complex_t r_i1 = c_i * ws.result_i(m), r_i2 = c_i * ws.result_i(m, true);
            return (r_o1 - r_i1) * F_arc * G_arc + (r_o2 - r_i2 - k2c_o * r_o1 + k2c_i * r_i1) * F * G * ws.n_dot_np(m);
        };
        for (int i = 0; i < Qtest; ++i) {
            for (int j = 0; j < Qtrial; ++j) {
                // Filling the matrix entry
                interaction_matrix(i, j) =
                -ii * 0.25 * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
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
        // The number of Reference Shape Functions in space
        int Q = space.getQ();
        DiscontinuousSpace<0> discont_space;
        double sqrtc_i = sqrt(c_i), sqrtc_o = sqrt(c_o);
        complex_t ksqrtc_i = k * sqrtc_i, ksqrtc_o = k * sqrtc_o, k2c_i = k * k * c_i, k2c_o = k * k * c_o;
        double ksqrtca_i = std::abs(ksqrtc_i), ksqrtca_o = std::abs(ksqrtc_o);
        unsigned N2 = N * N;
        bool with_i = c_i != 0.;
        // compute the required values once
        std::for_each (std::execution::par_unseq, ws.ind().cbegin(), ws.ind().cbegin() + N2, [&](auto &i) {
            unsigned I = i / N, J = i - I * N;
            ws.IJ(i) = std::make_pair(I, J);
            complex_t &result_i = ws.result_i(i), &result_o = ws.result_o(i), &result_ir = ws.result_i(i, true), &result_or = ws.result_o(i, true);
            result_i = result_o = result_ir = result_or = czero;
            double &t = ws.t(i), &s = ws.s(i), &w = ws.w(i), &t_norm = ws.t_norm(i), &tp_norm = ws.tp_norm(i), &n_dot_np = ws.n_dot_np(i);
            s = GaussQR.x(I);
            t = GaussQR.x(J);
            w = GaussQR.w(I) * GaussQR.w(J);
            // Finding the tangent of pi_p to get its normal
            const Eigen::Vector2d &tangent_p = pi_p.Derivative_01(t);
            const Eigen::Vector2d &tangent = pi.Derivative_01(s);
            tp_norm = tangent_p.norm();
            t_norm = tangent.norm();
            n_dot_np = (tangent_p(1) * tangent(1) + tangent_p(0) * tangent(0)) / (tp_norm * t_norm);
            double d = (pi[s] - pi_p[t]).norm();
            if (with_i && ksqrtca_i * d > epsilon ) {
                complex_t h0 = complex_bessel::H1(0, ksqrtc_i * d), h1 = complex_bessel::H1(1, ksqrtc_i * d);
                result_i = h0 * d * d - h1 * d / ksqrtc_i;
                result_ir = -h1 * d * 4.0 * ksqrtc_i + h0 * 2.0;
            }
            if (ksqrtca_o * d > epsilon ) {
                complex_t h0 = complex_bessel::H1(0, ksqrtc_o * d), h1 = complex_bessel::H1(1, ksqrtc_o * d);
                result_o = h0 * d * d - h1 * d / ksqrtc_o;
                result_or = -h1 * d * 4.0 * ksqrtc_o + h0 * 2.0;
            }
        });
        for (int i = 0; i < Q; ++i) {
            for (int j = 0; j < Q; ++j) {
                // Filling the matrix entry
                interaction_matrix(i, j) =
                -ii * 0.25 * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                    int m = ij.first * N + ij.second;
                    const double &s = ws.s(m), &t = ws.t(m);
                    auto F = space.evaluateShapeFunction(j, t) * ws.tp_norm(m);
                    auto G = space.evaluateShapeFunction(i, s) * ws.t_norm(m);
                    auto F_arc = space.evaluateShapeFunctionDot_01(j, t);
                    auto G_arc = space.evaluateShapeFunctionDot_01(i, s);
                    complex_t r_o1 = c_o * ws.result_o(m), r_o2 = c_o * ws.result_o(m, true);
                    complex_t r_i1 = c_i * ws.result_i(m), r_i2 = c_i * ws.result_i(m, true);
                    return sum + ws.w(m) * ((r_o1 - r_i1) * F_arc * G_arc + (r_o2 - r_i2 - k2c_o * r_o1 + k2c_i * r_i1) * F * G * ws.n_dot_np(m));
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
        else {// Disjoint panels case*/
            ComputeIntegralGeneral(interaction_matrix, pi, pi_p, space, GaussQR, k, c_i, c_o, ws);
        }
    }

    Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh &mesh,
                                    const AbstractBEMSpace &space,
                                    const QuadRule &GaussQR,
                                    const QuadRule &CGaussQR,
                                    const complex_t &k,
                                    double c_i,
                                    double c_o){
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
        Eigen::MatrixXcd output = Eigen::MatrixXd::Zero(dims, dims);
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

}// namespace hypersingular_helmholtz
