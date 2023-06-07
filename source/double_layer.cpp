#include "double_layer.hpp"
#include "cbessel.hpp"
#include <numeric>
#include <execution>

    namespace double_layer_helmholtz {

        typedef std::complex<double> complex_t;
        static const complex_t czero(0., 0.);
        static const double epsilon = std::numeric_limits<double>::epsilon();

        inline
        void ComputeIntegralCoinciding(Eigen::MatrixXcd &interaction_matrix,
                                       const AbstractParametrizedCurve &pi,
                                       const AbstractParametrizedCurve &pi_p,
                                       const AbstractBEMSpace &trial_space,
                                       const AbstractBEMSpace &test_space,
                                       const QuadRule &GaussQR,
                                       const complex_t &k,
                                       const double c_i, const double c_o,
                                       gq_workspace_t &ws) throw() {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = test_space.getQ();
            complex_t ksqrtc_i = k * sqrt(c_i), ksqrtc_o = k * sqrt(c_o);
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
                // Finding the tangent of pi_p to get its normal
                const Eigen::Vector2d &tangent_p = pi_p.Derivative_01(t);
                const Eigen::Vector2d &tangent = pi.Derivative_01(s);
                tp_norm = tangent_p.norm();
                t_norm = tangent.norm();
                Eigen::Vector2d normal;
                // Outward normal vector
                normal << tangent_p(1), -tangent_p(0);
                // Normalizing the normal vector
                normal = normal / tp_norm;
                Eigen::Vector2d v = pi[s] - pi_p[t];
                double v_norm = v.norm(), d = v.dot(normal) * 0.25 / v_norm;
                if (with_i) {
                    if (ksqrtca_i * v_norm > epsilon ) // Away from singularity
                        result_i = (ksqrtc_i * d) * complex_bessel::H1_1_i(ksqrtc_i * v_norm);
                    else if (v_norm > epsilon)
                        result_i = M_2_PI * d / v_norm;
                }
                if (ksqrtca_o * v_norm > epsilon ) // Away from singularity
                    result_o = (ksqrtc_o * d) * complex_bessel::H1_1_i(ksqrtc_o * v_norm);
                else if (v_norm > epsilon)
                    result_o = M_2_PI * d / v_norm;
            });
            // Lambda expression for the integrand
            auto integrand = [&](int i, int j, int m) {
                auto F = trial_space.evaluateShapeFunction(j, ws.t(m)) * ws.tp_norm(m);
                auto G = test_space.evaluateShapeFunction(i, ws.s(m)) * ws.t_norm(m);
                return (ws.result_o(m) - ws.result_i(m)) * F * G;
            };
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    interaction_matrix(i, j) =
                    std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                        int I = ij.first, J = ij.second;
                        return sum + ws.w(I * N + J) * (integrand(i, j, I * N + J) + integrand(i, j, J * N + I));
                    });
                }
            }
        }

        inline
        void ComputeIntegralAdjacent(Eigen::MatrixXcd &interaction_matrix,
                                     const AbstractParametrizedCurve &pi,
                                     const AbstractParametrizedCurve &pi_p,
                                     const AbstractBEMSpace &trial_space,
                                     const AbstractBEMSpace &test_space,
                                     const QuadRule &GaussQR,
                                     const complex_t &k,
                                     const double c_i, const double c_o,
                                     gq_workspace_t &ws) throw() {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = test_space.getQ();
            // Computing the (i,j)th matrix entry
            bool swap = (pi(1) - pi_p(-1)).norm() / 100. > epsilon;
            complex_t ksqrtc_i = k * sqrt(c_i), ksqrtc_o = k * sqrt(c_o);
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
                // Finding the tangent of pi_p to get its normal
                const Eigen::Vector2d &tangent_p = swap ? -pi_p.Derivative_01_swapped(t) : pi_p.Derivative_01(t);
                const Eigen::Vector2d &tangent = swap ? pi.Derivative_01(s) : pi.Derivative_01_swapped(s);
                tp_norm = tangent_p.norm();
                t_norm = tangent.norm();
                Eigen::Vector2d normal;
                normal << tangent_p(1), -tangent_p(0);
                // Normalizing the normal vector
                normal = normal / tp_norm;
                if (swap) {
                    Eigen::Vector2d v = pi[s] - pi_p.swapped_op(t);
                    double v_norm = v.norm(), d = v.dot(normal) * 0.25 / v_norm;
                    if (with_i) {
                        if (ksqrtca_i * v_norm > epsilon ) // Away from singularity
                            result_i = (ksqrtc_i * d) * complex_bessel::H1_1_i(ksqrtc_i * v_norm);
                        else if (v_norm > epsilon)
                            result_i = M_2_PI * d / v_norm;
                    }
                    if (ksqrtca_o * v_norm > epsilon ) // Away from singularity
                        result_o = (ksqrtc_o * d) * complex_bessel::H1_1_i(ksqrtc_o * v_norm);
                    else if (v_norm > epsilon)
                        result_o = M_2_PI * d / v_norm;
                } else {
                    Eigen::Vector2d v = pi.swapped_op(s) - pi_p[t];
                    double v_norm = v.norm(), d = v.dot(normal) * 0.25 / v_norm;
                    if (with_i) {
                        if (ksqrtca_i * v_norm > epsilon ) // Away from singularity
                            result_i = (ksqrtc_i * d) * complex_bessel::H1_1_i(ksqrtc_i * v_norm);
                        else if (v_norm > epsilon)
                            result_i = M_2_PI * d / v_norm;
                    }
                    if (ksqrtca_o * v_norm > epsilon ) // Away from singularity
                        result_o = (ksqrtc_o * d) * complex_bessel::H1_1_i(ksqrtc_o * v_norm);
                    else if (v_norm > epsilon)
                        result_o = M_2_PI * d / v_norm;
                }
            });
            // Lambda expression for the integrand
            auto integrand = [&](int i, int j, int m) {
                const double &s = ws.s(m), &t = ws.t(m);
                auto F = (swap ? trial_space.evaluateShapeFunction_01_swapped(j, t) : trial_space.evaluateShapeFunction(j, t)) * ws.tp_norm(m);
                auto G = (swap ? test_space.evaluateShapeFunction(i, s) : test_space.evaluateShapeFunction_01_swapped(i, s)) * ws.t_norm(m);
                return (ws.result_o(m) - ws.result_i(m)) * F * G;
            };
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    interaction_matrix(i, j) =
                    std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                        int I = ij.first, J = ij.second;
                        return sum + ws.w(I * N + J) * (integrand(i, j, I * N + J) + integrand(i, j, J * N + I));
                    });
                }
            }
        }

        inline
        void ComputeIntegralGeneral(Eigen::MatrixXcd &interaction_matrix,
                                    const AbstractParametrizedCurve &pi,
                                    const AbstractParametrizedCurve &pi_p,
                                    const AbstractBEMSpace &trial_space,
                                    const AbstractBEMSpace &test_space,
                                    const QuadRule &GaussQR,
                                    const complex_t &k,
                                    const double c_i, const double c_o,
                                    gq_workspace_t &ws) throw() {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in space
            int Qtest = test_space.getQ();
            // Computing the (i,j)th matrix entry
            complex_t ksqrtc_i = k * sqrt(c_i), ksqrtc_o = k * sqrt(c_o);
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
                // Finding the tangent of pi_p to get its normal
                const Eigen::Vector2d &tangent = pi_p.Derivative_01(t);
                const Eigen::Vector2d &tangent_p = pi.Derivative_01(s);
                t_norm = tangent.norm();
                tp_norm = tangent_p.norm();
                Eigen::Vector2d normal;
                // Outward normal vector
                normal << tangent(1), -tangent(0);
                // Normalizing the normal vector
                normal = normal / t_norm;
                Eigen::Vector2d v = pi[s] - pi_p[t];
                double v_norm = v.norm(), d = v.dot(normal) * 0.25 / v_norm;
                if (with_i) {
                    if (ksqrtca_i * v_norm > epsilon ) // Away from singularity
                        result_i = (ksqrtc_i * d) * complex_bessel::H1_1_i(ksqrtc_i * v_norm);
                    else if (v_norm > epsilon)
                        result_i = M_2_PI * d / v_norm;
                }
                if (ksqrtca_o * v_norm > epsilon ) // Away from singularity
                    result_o = (ksqrtc_o * d) * complex_bessel::H1_1_i(ksqrtc_o * v_norm);
                else if (v_norm > epsilon)
                    result_o = M_2_PI * d / v_norm;
            });
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    interaction_matrix(i, j) =
                    std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                        int m = ij.first * N + ij.second;
                        auto F = trial_space.evaluateShapeFunction(j, ws.t(m)) * ws.t_norm(m);
                        auto G = test_space.evaluateShapeFunction(i, ws.s(m)) * ws.tp_norm(m);
                        return sum + ws.w(m) * (ws.result_o(m) - ws.result_i(m)) * F * G;
                    });
                }
            }
        }

        inline
        void InteractionMatrix(Eigen::MatrixXcd &interaction_matrix,
                               const AbstractParametrizedCurve &pi,
                               const AbstractParametrizedCurve &pi_p,
                               const AbstractBEMSpace &trial_space,
                               const AbstractBEMSpace &test_space,
                               const QuadRule &GaussQR,
                               const QuadRule &CGaussQR,
                               const complex_t &k,
                               const double c_i, const double c_o,
                               gq_workspace_t &ws) throw() {
            if (&pi == &pi_p) { // Same Panels case
                ComputeIntegralCoinciding(interaction_matrix, pi, pi_p, trial_space, test_space, CGaussQR, k, c_i, c_o, ws);
            }
            else if ((pi(1) - pi_p(-1)).norm() / 100. < epsilon ||
                     (pi(-1) - pi_p(1)).norm() / 100. < epsilon) {// Adjacent Panels case
                ComputeIntegralAdjacent(interaction_matrix, pi, pi_p, trial_space, test_space, CGaussQR, k, c_i, c_o, ws);
            }
            else { //Disjoint panels case
                ComputeIntegralGeneral(interaction_matrix, pi, pi_p, trial_space, test_space, GaussQR, k, c_i, c_o, ws);
            }
        }

        Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                        const AbstractBEMSpace &trial_space,
                                        const AbstractBEMSpace &test_space,
                                        const unsigned int &N,
                                        const complex_t &k,
                                        const double c_i,
                                        const double c_o) {
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
            Eigen::MatrixXcd interaction_matrix(Qtest, Qtrial);
            // Initializing the Galerkin matrix with zeros
            Eigen::MatrixXcd output = Eigen::MatrixXcd::Zero(rows, cols);
            // Panel oriented assembly
            QuadRule GaussQR = getGaussQR(N,0.,1.);
            QuadRule CGaussQR = getCGaussQR(N);
            gq_workspace_t ws(std::max(GaussQR.n, CGaussQR.n));
            unsigned int i, j, I, J;
            for (i = 0; i < numpanels; ++i) {
                for (j = 0; j < numpanels; ++j) {
                    // Getting the interaction matrix for the pair of panels i and j
                    InteractionMatrix(interaction_matrix, *panels[i], *panels[j], trial_space, test_space, GaussQR, CGaussQR, k, c_i, c_o, ws);
                    // Local to global mapping of the elements in interaction matrix
                    for (I = 0; I < Qtest; ++I) {
                        for (J = 0; J < Qtrial; ++J) {
                            // Filling the Galerkin matrix entries
                            output(test_space.LocGlobMap(I + 1, i + 1, numpanels) - 1, trial_space.LocGlobMap(J + 1, j + 1, numpanels) - 1)
                                += interaction_matrix(I, J);
                        }
                    }
                }
            }
            return output;
        }

    } // namespace double_layer_helmholtz
