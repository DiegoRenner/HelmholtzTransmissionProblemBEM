#include "double_layer_der.hpp"
#include "cbessel.h"
#include <execution>

    namespace double_layer_helmholtz_der {

        typedef std::complex<double> complex_t;
        static const complex_t ii(0., 1.);
        static const complex_t czero(0., 0.);
        static const double epsilon = std::numeric_limits<double>::epsilon();

        void InteractionMatrix(Eigen::MatrixXcd &interaction_matrix,
                               const AbstractParametrizedCurve &pi,
                               const AbstractParametrizedCurve &pi_p,
                               const AbstractBEMSpace &trial_space,
                               const AbstractBEMSpace &test_space,
                               const QuadRule &GaussQR,
                               const QuadRule &CGaussQR,
                               const complex_t &k,
                               const double c_i, const double c_o,
                               gq_workspace_t &ws) {
            if (&pi == &pi_p) { // Same Panels case
                ComputeIntegralCoinciding(interaction_matrix, pi, pi_p, trial_space, test_space, CGaussQR, k, c_i, c_o, ws);
            }
            else if ((pi(1) - pi_p(-1)).norm() / 100. < epsilon ||
                     (pi(-1) - pi_p(1)).norm() / 100. < epsilon) {// Adjacent Panels case
                ComputeIntegralAdjacent(interaction_matrix, pi, pi_p, trial_space, test_space, CGaussQR, k, c_i, c_o, ws);
            } else { //Disjoint panels case
                ComputeIntegralGeneral(interaction_matrix, pi, pi_p, trial_space, test_space, GaussQR, k, c_i, c_o, ws);
            }
        }

        void ComputeIntegralCoinciding(Eigen::MatrixXcd &interaction_matrix,
                                       const AbstractParametrizedCurve &pi,
                                       const AbstractParametrizedCurve &pi_p,
                                       const AbstractBEMSpace &trial_space,
                                       const AbstractBEMSpace &test_space,
                                       const QuadRule &GaussQR,
                                       const complex_t &k,
                                       const double c_i, const double c_o,
                                       gq_workspace_t &ws) {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = test_space.getQ();
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
                Eigen::Vector2d v = pi[s]-pi_p[t];
                double d = v.norm(), vdotn = v.dot(normal);
                if (with_i && ksqrtca_i * d > epsilon ) // Away from singularity
                    result_i = ComplexBessel::H1(0, ksqrtc_i * d) * ksqrtc_i * vdotn;
                if (ksqrtca_o * d > epsilon ) // Away from singularity
                    result_o = ComplexBessel::H1(0, ksqrtc_o * d) * ksqrtc_o * vdotn;
            });
            // Lambda expression for the integrand
            auto integrand = [&](int i, int j, int m) {
                auto F = trial_space.evaluateShapeFunction(j, ws.t(m)) * ws.tp_norm(m);
                auto G = test_space.evaluateShapeFunction(i, ws.s(m)) * ws.t_norm(m);
                return (sqrtc_o * ws.result_o(m) - sqrtc_i * ws.result_i(m)) * F * G;
            };
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    interaction_matrix(i, j) =
                    0.25 * ii * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                        const int &I = ij.first, &J = ij.second;
                        return sum + ws.w(I * N + J) * (integrand(i, j, I * N + J) + integrand(i, j, J * N + I));
                    });
                }
            }
        }

        void ComputeIntegralAdjacent(Eigen::MatrixXcd &interaction_matrix,
                                     const AbstractParametrizedCurve &pi,
                                     const AbstractParametrizedCurve &pi_p,
                                     const AbstractBEMSpace &trial_space,
                                     const AbstractBEMSpace &test_space,
                                     const QuadRule &GaussQR,
                                     const complex_t &k,
                                     const double c_i, const double c_o,
                                     gq_workspace_t &ws) {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in trial space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in test space
            int Qtest = test_space.getQ();
            bool swap = ((pi(1) - pi_p(-1)).norm() / 100. > sqrt(epsilon));
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
                // Finding the tangent of pi_p to get its normal
                const Eigen::Vector2d &tangent_p = swap ? pi_p.Derivative_01_swapped(t) : pi_p.Derivative_01(t);
                const Eigen::Vector2d &tangent = swap ? pi.Derivative_01(s) : pi.Derivative_01_swapped(s);
                tp_norm = tangent_p.norm();
                t_norm = tangent.norm();
                Eigen::Vector2d normal;
                // Outward normal vector
                if (swap) {
                    normal << -tangent_p(1), tangent_p(0);
                } else {
                    normal << tangent_p(1), -tangent_p(0);
                }
                // Normalizing the normal vector
                normal = normal / tp_norm;
                Eigen::Vector2d v = swap ? pi[s]-pi_p.swapped_op(t) : pi.swapped_op(s)-pi_p[t];
                double d = v.norm(), vdotn = v.dot(normal);
                if (with_i && ksqrtca_i * d > epsilon)
                    result_i = ComplexBessel::H1(0, ksqrtc_i * d) * ksqrtc_i * vdotn;
                if (ksqrtca_o * d > epsilon)
                    result_o = ComplexBessel::H1(0, ksqrtc_o * d) * ksqrtc_o * vdotn;
            });
            // Lambda expression for the integrand
            auto integrand = [&](int i, int j, int m) {
                auto F = (swap ? trial_space.evaluateShapeFunction_01_swapped(j, ws.t(m)) : trial_space.evaluateShapeFunction(j, ws.t(m))) * ws.tp_norm(m);
                auto G = (swap ? test_space.evaluateShapeFunction(i, ws.s(m)) : test_space.evaluateShapeFunction_01_swapped(i, ws.s(m))) * ws.t_norm(m);
                return (sqrtc_o * ws.result_o(m) - sqrtc_i * ws.result_i(m)) * F * G;
            };
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    interaction_matrix(i, j) =
                    0.25 * ii * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                        const int &I = ij.first, &J = ij.second;
                        return sum + ws.w(I * N + J) * (integrand(i, j, I * N + J) + integrand(i, j, J * N + I));
                    });
                }
            }
        }

        void ComputeIntegralGeneral(Eigen::MatrixXcd &interaction_matrix,
                                    const AbstractParametrizedCurve &pi,
                                    const AbstractParametrizedCurve &pi_p,
                                    const AbstractBEMSpace &trial_space,
                                    const AbstractBEMSpace &test_space,
                                    const QuadRule &GaussQR,
                                    const complex_t &k,
                                    const double c_i, const double c_o,
                                    gq_workspace_t &ws) {
            unsigned N = GaussQR.n; // quadrature order for the GaussQR object.
            // The number of Reference Shape Functions in space
            int Qtrial = trial_space.getQ();
            // The number of Reference Shape Functions in space
            int Qtest = test_space.getQ();
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
                Eigen::Vector2d v = pi[s]-pi_p[t];
                double d = v.norm(), vdotn = v.dot(normal);
                if (with_i && ksqrtca_i * d > epsilon)
                    result_i = ComplexBessel::H1(0, ksqrtc_i * d) * ksqrtc_i * vdotn;
                if (ksqrtca_o * d > epsilon)
                    result_o = ComplexBessel::H1(0, ksqrtc_o * d) * ksqrtc_o * vdotn;
            });
            for (int i = 0; i < Qtest; ++i) {
                for (int j = 0; j < Qtrial; ++j) {
                    interaction_matrix(i, j) =
                    0.25 * ii * std::accumulate(ws.IJ().cbegin(), ws.IJ().cbegin() + N2, czero, [&](const auto &sum, const auto &ij) {
                        int m = ij.first * N + ij.second;
                        auto F = trial_space.evaluateShapeFunction(j, ws.t(m)) * ws.tp_norm(m);
                        auto G = test_space.evaluateShapeFunction(i, ws.s(m)) * ws.t_norm(m);
                        return sum + ws.w(m) * (sqrtc_o * ws.result_o(m) - sqrtc_i * ws.result_i(m)) * F * G;
                    });
                }
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
            Eigen::MatrixXcd output = Eigen::MatrixXd::Zero(rows, cols);
            // Panel oriented assembly
            QuadRule GaussQR = getGaussQR(N,0.,1.);
            QuadRule CGaussQR = getCGaussQR(N);
            gq_workspace_t ws(std::max(GaussQR.n, CGaussQR.n));
            unsigned i, j, I, J;
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

