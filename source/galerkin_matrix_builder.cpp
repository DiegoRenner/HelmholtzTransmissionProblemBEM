#include "galerkin_matrix_builder.hpp"
#include "cbessel.hpp"
#include <numeric>
#include <iostream>
#include <chrono>

using namespace std::complex_literals;

typedef std::complex<double> complex_t;
static const complex_t czero(0., 0.);
static const double epsilon = std::numeric_limits<double>::epsilon();

GalerkinMatrixBuilder::GalerkinMatrixBuilder(const ParametrizedMesh &mesh,
                                             const AbstractBEMSpace &test_space_in,
                                             const AbstractBEMSpace &trial_space_in,
                                             const QuadRule &GaussQR_in,
                                             const QuadRule &CGaussQR_in)
: test_space(test_space_in), trial_space(trial_space_in), panels(mesh.getPanels()), GaussQR(GaussQR_in), CGaussQR(CGaussQR_in)
{
    hankel_time = mapping_time = interaction_matrix_time = 0;
    numpanels = mesh.getNumPanels();
    dim_test = test_space.getSpaceDim(numpanels);
    dim_trial = trial_space.getSpaceDim(numpanels);
    Qtest = test_space.getQ();
    Qtrial = trial_space.getQ();
    size_t N = CGaussQR.n, Ns = GaussQR.n;
    m_h0.resize(N, 2 * N);
    m_h1.resize(N, 2 * N);
    m_h_arg.resize(N, 2 * N);
    m_h_arg_half.resize(N, N);
    m_h_arg_small.resize(Ns, Ns);
    m_h0_res_half.resize(N, N);
    m_h1_res_half.resize(N, N);
    m_h0_res_small.resize(Ns, Ns);
    m_h1_res_small.resize(Ns, Ns);
    m_v.resize(N, 2 * N);
    m_v_norm.resize(N, 2 * N);
    m_tangent.resize(N, 2 * N);
    m_tangent_p.resize(N, 2 * N);
    m_tangent_norm.resize(N, 2 * N);
    m_tangent_p_norm.resize(N, 2 * N);
    m_zero.setZero(N, N);
    m_cf.resize(N, N);
    m_vdotn.resize(N, N);
    m_temp.resize(N, N);
    m_tc.resize(N, N);
    m_ta.resize(N, N);
    m_tg.resize(Ns, Ns);
    m_sc.resize(N, N);
    m_sa.resize(N, N);
    m_sg.resize(Ns, Ns);
    m_wc.resize(N, N);
    m_wa.resize(N, N);
    m_wg.resize(Ns, Ns);
    for (size_t I = 0; I < N; ++I) {
        for (size_t J = 0; J < N; ++J) {
            double t = CGaussQR.x(J), s = t * (1. - CGaussQR.x(I)), w = t * CGaussQR.w(I) * CGaussQR.w(J);
            m_tc(I, J) = t; m_sc(I, J) = s; m_wc(I, J) = w;
            t = CGaussQR.x(I), s = t * CGaussQR.x(J), w = t * CGaussQR.w(I) * CGaussQR.w(J);
            m_ta(I, J) = t; m_sa(I, J) = s; m_wa(I, J) = w;
            if (I < Ns && J < Ns) {
                t = GaussQR.x(J), s = GaussQR.x(I), w = GaussQR.w(I) * GaussQR.w(J);
                m_tg(I, J) = t; m_sg(I, J) = s; m_wg(I, J) = w;
            }
        }
    }
    m_double_layer_coinciding_fg.resize(Qtest * N, Qtrial * N);
    m_double_layer_coinciding_fg_t.resize(Qtest * N, Qtrial * N);
    m_double_layer_adjacent_fg.resize(Qtest * N, Qtrial * N);
    m_double_layer_adjacent_fg_swap.resize(Qtest * N, Qtrial * N);
    m_double_layer_adjacent_fg_t.resize(Qtest * N, Qtrial * N);
    m_double_layer_adjacent_fg_swap_t.resize(Qtest * N, Qtrial * N);
    m_double_layer_general_fg.resize(Qtest * Ns, Qtrial * Ns);
    m_double_layer_general_fg_t.resize(Qtest * Ns, Qtrial * Ns);
    m_hypersingular_coinciding_fg.resize(Qtrial * N, Qtrial * N);
    m_hypersingular_coinciding_fg_arc.resize(Qtrial * N, Qtrial * N);
    m_hypersingular_adjacent_fg.resize(Qtrial * N, Qtrial * N);
    m_hypersingular_adjacent_fg_swap.resize(Qtrial * N, Qtrial * N);
    m_hypersingular_adjacent_fg_arc.resize(Qtrial * N, Qtrial * N);
    m_hypersingular_adjacent_fg_arc_swap.resize(Qtrial * N, Qtrial * N);
    m_hypersingular_general_fg.resize(Qtrial * Ns, Qtrial * Ns);
    m_hypersingular_general_fg_arc.resize(Qtrial * Ns, Qtrial * Ns);
    m_single_layer_coinciding_fg.resize(Qtest * N, Qtest * N);
    m_single_layer_adjacent_fg.resize(Qtest * N, Qtest * N);
    m_single_layer_adjacent_fg_swap.resize(Qtest * N, Qtest * N);
    m_single_layer_general_fg.resize(Qtest * Ns, Qtest * Ns);
    for (size_t j = 0; j < Qtrial; ++j) {
        for (size_t i = 0; i < Qtest; ++i) {
            auto dbl_g_fg = m_double_layer_general_fg.block(i * Ns, j * Ns, Ns, Ns);
            auto dbl_g_fg_t = m_double_layer_general_fg_t.block(i * Ns, j * Ns, Ns, Ns);
            for (size_t J = 0; J < N; ++J) {
                for (size_t I = 0; I < N; ++I) {
                    double sc = m_sc(I, J), sa = m_sa(I, J), tc = m_tc(I, J), ta = m_ta(I, J);
                    m_double_layer_coinciding_fg(i * N + I, j * N + J) =
                        test_space.evaluateShapeFunction(i, sc) * trial_space.evaluateShapeFunction(j, tc);
                    m_double_layer_coinciding_fg_t(i * N + I, j * N + J) =
                        test_space.evaluateShapeFunction(i, tc) * trial_space.evaluateShapeFunction(j, sc);
                    m_double_layer_adjacent_fg(i * N + I, j * N + J) =
                        test_space.evaluateShapeFunction_01_swapped(i, sa) * trial_space.evaluateShapeFunction(j, ta);
                    m_double_layer_adjacent_fg_swap(i * N + I, j * N + J) =
                        test_space.evaluateShapeFunction(i, sa) * trial_space.evaluateShapeFunction_01_swapped(j, ta);
                    m_double_layer_adjacent_fg_t(i * N + I, j * N + J) =
                        test_space.evaluateShapeFunction_01_swapped(i, ta) * trial_space.evaluateShapeFunction(j, sa);
                    m_double_layer_adjacent_fg_swap_t(i * N + I, j * N + J) =
                        test_space.evaluateShapeFunction(i, ta) * trial_space.evaluateShapeFunction_01_swapped(j, sa);
                    if (I < Ns && J < Ns) {
                        double s = m_sg(I, J), t = m_tg(I, J);
                        dbl_g_fg(I, J) =
                            test_space.evaluateShapeFunction(i, s) * trial_space.evaluateShapeFunction(j, t);
                        dbl_g_fg_t(I, J) =
                            test_space.evaluateShapeFunction(i, t) * trial_space.evaluateShapeFunction(j, s);
                    }
                }
            }
        }
        for (size_t i = 0; i < Qtrial; ++i) {
            auto hyp_c_fg = m_hypersingular_coinciding_fg.block(i * N, j * N, N, N);
            auto hyp_c_fg_arc = m_hypersingular_coinciding_fg_arc.block(i * N, j * N, N, N);
            auto hyp_a_fg = m_hypersingular_adjacent_fg.block(i * N, j * N, N, N);
            auto hyp_a_fg_swap = m_hypersingular_adjacent_fg_swap.block(i * N, j * N, N, N);
            auto hyp_a_fg_arc = m_hypersingular_adjacent_fg_arc.block(i * N, j * N, N, N);
            auto hyp_a_fg_arc_swap = m_hypersingular_adjacent_fg_arc_swap.block(i * N, j * N, N, N);
            auto hyp_g_fg = m_hypersingular_general_fg.block(i * Ns, j * Ns, Ns, Ns);
            auto hyp_g_fg_arc = m_hypersingular_general_fg_arc.block(i * Ns, j * Ns, Ns, Ns);
            for (size_t J = 0; J < N; ++J) {
                for (size_t I = 0; I < N; ++I) {
                    double sc = m_sc(I, J), sa = m_sa(I, J), tc = m_tc(I, J), ta = m_ta(I, J);
                    hyp_c_fg(I, J) =
                        trial_space.evaluateShapeFunction(i, sc) * trial_space.evaluateShapeFunction(j, tc);
                    hyp_c_fg_arc(I, J) =
                        trial_space.evaluateShapeFunctionDot_01(i, sc) * trial_space.evaluateShapeFunctionDot_01(j, tc);
                    hyp_a_fg(I, J) =
                        trial_space.evaluateShapeFunction_01_swapped(i, sa) * trial_space.evaluateShapeFunction(j, ta);
                    hyp_a_fg_swap(I, J) =
                        trial_space.evaluateShapeFunction(i, sa) * trial_space.evaluateShapeFunction_01_swapped(j, ta);
                    hyp_a_fg_arc(I, J) =
                        trial_space.evaluateShapeFunctionDot_01_swapped(i, sa) * trial_space.evaluateShapeFunctionDot_01(j, ta);
                    hyp_a_fg_arc_swap(I, J) =
                        trial_space.evaluateShapeFunctionDot_01(i, sa) * trial_space.evaluateShapeFunctionDot_01_swapped(j, ta);
                    if (I < Ns && J < Ns) {
                        double s = m_sg(I, J), t = m_tg(I, J);
                        hyp_g_fg(I, J) =
                            trial_space.evaluateShapeFunction(i, s) * trial_space.evaluateShapeFunction(j, t);
                        hyp_g_fg_arc(I, J) =
                            trial_space.evaluateShapeFunctionDot_01(i, s) * trial_space.evaluateShapeFunctionDot_01(j, t);
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < Qtest; ++i) {
        for (size_t j = 0; j < Qtest; ++j) {
            auto sng_c_fg = m_single_layer_coinciding_fg.block(i * N, j * N, N, N);
            auto sng_a_fg = m_single_layer_adjacent_fg.block(i * N, j * N, N, N);
            auto sng_a_fg_swap = m_single_layer_adjacent_fg_swap.block(i * N, j * N, N, N);
            auto sng_g_fg = m_single_layer_general_fg.block(i * Ns, j * Ns, Ns, Ns);
            for (size_t J = 0; J < N; ++J) {
                for (size_t I = 0; I < N; ++I) {
                    double sc = m_sc(I, J), sa = m_sa(I, J), tc = m_tc(I, J), ta = m_ta(I, J);
                    sng_c_fg(I, J) =
                        test_space.evaluateShapeFunction(i, sc) * test_space.evaluateShapeFunction(j, tc);
                    sng_a_fg(I, J) =
                        test_space.evaluateShapeFunction_01_swapped(i, sa) * test_space.evaluateShapeFunction(j, ta);
                    sng_a_fg_swap(I, J) =
                        test_space.evaluateShapeFunction(i, sa) * test_space.evaluateShapeFunction_01_swapped(j, ta);
                    if (I < Ns && J < Ns) {
                        double s = m_sg(I, J), t = m_tg(I, J);
                        sng_g_fg(I, J) =
                            test_space.evaluateShapeFunction(i, s) * test_space.evaluateShapeFunction(j, t);
                    }
                }
            }
        }
    }
    indN2.resize(2 * N);
    std::iota(indN2.begin(), indN2.end(), 0);
    indNs.resize(Ns);
    std::iota(indNs.begin(), indNs.end(), 0);
    double_layer_interaction_matrix.resize(Qtest, Qtrial);
    double_layer_der_interaction_matrix.resize(Qtest, Qtrial);
    double_layer_der2_interaction_matrix.resize(Qtest, Qtrial);
    hypersingular_interaction_matrix.resize(Qtrial, Qtrial);
    hypersingular_der_interaction_matrix.resize(Qtrial, Qtrial);
    hypersingular_der2_interaction_matrix.resize(Qtrial, Qtrial);
    single_layer_interaction_matrix.resize(Qtest, Qtest);
    single_layer_der_interaction_matrix.resize(Qtest, Qtest);
    single_layer_der2_interaction_matrix.resize(Qtest, Qtest);
    double_layer_matrix.resize(dim_test, dim_trial);
    double_layer_der_matrix.resize(dim_test, dim_trial);
    double_layer_der2_matrix.resize(dim_test, dim_trial);
    hypersingular_matrix.resize(dim_trial, dim_trial);
    hypersingular_der_matrix.resize(dim_trial, dim_trial);
    hypersingular_der2_matrix.resize(dim_trial, dim_trial);
    single_layer_matrix.resize(dim_test, dim_test);
    single_layer_der_matrix.resize(dim_test, dim_test);
    single_layer_der2_matrix.resize(dim_test, dim_test);
}

// compute values required for the coinciding panels case
void GalerkinMatrixBuilder::compute_coinciding(const AbstractParametrizedCurve &p) throw() {
    size_t N = CGaussQR.n;
    for (size_t J = 0; J < N; ++J) {
        for (size_t I = 0; I < N; ++I) {
            double t = m_tc(I, J), s = m_sc(I, J);
            Eigen::Vector2d tangent_p = p.Derivative_01(t);
            Eigen::Vector2d tangent = p.Derivative_01(s);
            Eigen::Vector2d v = p[s] - p[t];
            double v_norm = v.norm();
            complex_t arg = ksqrtc * v_norm;
            m_tangent_p(I, J) = complex_t(tangent_p(0), tangent_p(1));
            m_tangent(I, J) = complex_t(tangent(0), tangent(1));
            m_tangent_p_norm(I, J) = tangent_p.norm();
            m_tangent_norm(I, J) = tangent.norm();
            m_v(I, J) = complex_t(v(0), v(1));
            m_v_norm(I, J) = v_norm;
            if (!k_real_positive) {
                m_h0(I, J) = 1i * complex_bessel::H1(0, arg);
                m_h1(I, J) = 1i * complex_bessel::H1(1, arg);
            } else m_h_arg_half(I, J) = arg.real();
        }
    }
    if (k_real_positive) {
        auto tic = chrono::high_resolution_clock::now();
        complex_bessel::H1_01_i(m_h_arg_half, m_h0_res_half, m_h1_res_half);
        auto toc = chrono::high_resolution_clock::now();
        auto dur = chrono::duration_cast<std::chrono::microseconds>(toc - tic);
        hankel_time += dur.count();
        m_h0.block(0, 0, N, N) = m_h0_res_half;
        m_h1.block(0, 0, N, N) = m_h1_res_half;
    }
}

// compute values required for the adjacent panels case
void GalerkinMatrixBuilder::compute_adjacent(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p, bool swap) throw() {
    size_t N = CGaussQR.n;
    for (size_t J = 0; J < 2 * N; ++J) {
        for (size_t I = 0; I < N; ++I) {
            double t = m_ta(I, J % N), s = m_sa(I, J % N);
            if (J >= N)
                std::swap(s, t);
            Eigen::Vector2d tangent_p = swap ? -pi_p.Derivative_01_swapped(t) : pi_p.Derivative_01(t);
            Eigen::Vector2d tangent = swap ? pi.Derivative_01(s) : -pi.Derivative_01_swapped(s);
            Eigen::Vector2d v = swap ? pi[s] - pi_p.swapped_op(t) : pi.swapped_op(s) - pi_p[t];
            double v_norm = v.norm();
            complex_t arg = ksqrtc * v_norm;
            m_tangent_p(I, J) = complex_t(tangent_p(0), tangent_p(1));
            m_tangent(I, J) = complex_t(tangent(0), tangent(1));
            m_tangent_p_norm(I, J) = tangent_p.norm();
            m_tangent_norm(I, J) = tangent.norm();
            m_v(I, J) = complex_t(v(0), v(1));
            m_v_norm(I, J) = v_norm;
            if (!k_real_positive) {
                m_h0(I, J) = 1i * complex_bessel::H1(0, arg);
                m_h1(I, J) = 1i * complex_bessel::H1(1, arg);
            } else m_h_arg(I, J) = arg.real();
        }
    }
    if (k_real_positive) {
        auto tic = chrono::high_resolution_clock::now();
        complex_bessel::H1_01_i(m_h_arg, m_h0, m_h1);
        auto toc = chrono::high_resolution_clock::now();
        auto dur = chrono::duration_cast<std::chrono::microseconds>(toc - tic);
        hankel_time += dur.count();
    }
}

// compute values required for the disjoint panels case
void GalerkinMatrixBuilder::compute_general(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p) throw() {
    size_t N = GaussQR.n;
    auto tangent_p_block = m_tangent_p.block(0, 0, N, N);
    auto tangent_block = m_tangent.block(0, 0, N, N);
    auto tangent_p_norm_block = m_tangent_p_norm.block(0, 0, N, N);
    auto tangent_norm_block = m_tangent_norm.block(0, 0, N, N);
    auto v_block = m_v.block(0, 0, N, N);
    auto v_norm_block = m_v_norm.block(0, 0, N, N);
    auto h0_block = m_h0.block(0, 0, N, N);
    auto h1_block = m_h1.block(0, 0, N, N);
    for (size_t J = 0; J < N; ++J) {
        for (size_t I = 0; I < N; ++I) {
            double t = m_tg(I, J), s = m_sg(I, J);
            Eigen::Vector2d tangent_p = pi_p.Derivative_01(t);
            Eigen::Vector2d tangent = pi.Derivative_01(s);
            Eigen::Vector2d v = pi[s] - pi_p[t];
            double v_norm = v.norm();
            complex_t arg = ksqrtc * v_norm;
            tangent_p_block(I, J) = complex_t(tangent_p(0), tangent_p(1));
            tangent_block(I, J) = complex_t(tangent(0), tangent(1));
            tangent_p_norm_block(I, J) = tangent_p.norm();
            tangent_norm_block(I, J) = tangent.norm();
            v_block(I, J) = complex_t(v(0), v(1));
            v_norm_block(I, J) = v_norm;
            if (!k_real_positive) {
                h0_block(I, J) = 1i * complex_bessel::H1(0, arg);
                h1_block(I, J) = 1i * complex_bessel::H1(1, arg);
            } else m_h_arg_small(I, J) = arg.real();
        }
    }
    if (k_real_positive) {
        auto tic = chrono::high_resolution_clock::now();
        complex_bessel::H1_01_i(m_h_arg_small, m_h0_res_small, m_h1_res_small);
        auto toc = chrono::high_resolution_clock::now();
        auto dur = chrono::duration_cast<std::chrono::microseconds>(toc - tic);
        hankel_time += dur.count();
        h0_block = m_h0_res_small;
        h1_block = m_h1_res_small;
    }
}

// compute interaction matrices for the K submatrix, coinciding panels case
inline void GalerkinMatrixBuilder::double_layer_coinciding(int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    auto v = m_v.block(0, 0, N, N);
    auto v_norm = m_v_norm.block(0, 0, N, N);
    auto h1 = m_h1.block(0, 0, N, N), h0 = m_h0.block(0, 0, N, N);
    auto mask1 = (ksqrtca * v_norm) > epsilon;
    auto mask2 = v_norm > epsilon;
    auto otherwise = mask2.select(M_2_PI / (v_norm * v_norm), m_zero);
    m_temp = h0 - h1 * ksqrtc * v_norm;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (K = 0; K < 2; ++K) {
                auto fg = (K == 0 ? m_double_layer_coinciding_fg : m_double_layer_coinciding_fg_t).block(i * N, j * N, N, N);
                auto tangent = (K == 0 ? m_tangent_p : m_tangent).block(0, 0, N, N);
                auto t_norm = (K == 0 ? m_tangent_norm : m_tangent_p_norm).block(0, 0, N, N);
                m_vdotn = (K * 2 - 1) * (v.imag() * tangent.real() - v.real() * tangent.imag());
                m_cf = m_wc * t_norm * m_vdotn * fg;
                integral += (m_cf * mask1.select(ksqrtc * h1 / v_norm, otherwise)).sum();
                if (der > 0)
                    integral_der += ksqrtc * mask1.select(h0 * m_cf, m_zero).sum();
                if (der > 1)
                    integral_der2 += mask1.select(m_temp * m_cf, m_zero).sum();
            }
            double_layer_interaction_matrix(i, j) = integral;
            if (der > 0)
                double_layer_der_interaction_matrix(i, j) = sqrtc * integral_der;
            if (der > 1)
                double_layer_der2_interaction_matrix(i, j) = c * integral_der2;
        }
    }
}

// compute interaction matrices for the K submatrix, adjacent panels case
inline void GalerkinMatrixBuilder::double_layer_adjacent(bool swap, int der, bool transp) throw() {
    size_t i, j, K, N = CGaussQR.n;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (K = 0; K < 2; ++K) {
                auto fg = (swap ? ((transp ? K == 0 : K == 1) ? m_double_layer_adjacent_fg_swap_t
                                                              : m_double_layer_adjacent_fg_swap)
                                : ((transp ? K == 0 : K == 1) ? m_double_layer_adjacent_fg_t
                                                              : m_double_layer_adjacent_fg)).block(i * N, j * N, N, N);
                auto tangent = (transp ? m_tangent : m_tangent_p).block(0, K * N, N, N);
                auto t_norm = (transp ? m_tangent_p_norm : m_tangent_norm).block(0, K * N, N, N);
                auto v = m_v.block(0, K * N, N, N);
                auto v_norm = m_v_norm.block(0, K * N, N, N);
                auto h1 = m_h1.block(0, K * N, N, N), h0 = m_h0.block(0, K * N, N, N);
                auto mask = (ksqrtca * v_norm) > epsilon;
                m_vdotn = (transp? 1.0 : -1.0) * (v.imag() * tangent.real() - v.real() * tangent.imag());
                m_cf = m_wa * t_norm * m_vdotn * fg;
                integral += (m_cf * mask.select(ksqrtc * h1 / v_norm, (v_norm > epsilon).select(M_2_PI / (v_norm * v_norm), m_zero))).sum();
                if (der > 0)
                    integral_der += ksqrtc * mask.select(h0 * m_cf, m_zero).sum();
                if (der > 1)
                    integral_der2 += mask.select((h0 - h1 * ksqrtc * v_norm) * m_cf, m_zero).sum();
            }
            double_layer_interaction_matrix(i, j) = integral;
            if (der > 0)
                double_layer_der_interaction_matrix(i, j) = sqrtc * integral_der;
            if (der > 1)
                double_layer_der2_interaction_matrix(i, j) = c * integral_der2;
        }
    }
}

// compute interaction matrices for the K submatrix, disjoint panels case
inline void GalerkinMatrixBuilder::double_layer_general(int der, bool transp) throw() {
    size_t i, j, N = GaussQR.n;
    auto tangent = (transp ? m_tangent : m_tangent_p).block(0, 0, N, N);
    auto v = m_v.block(0, 0, N, N);
    auto t_norm = (transp ? m_tangent_p_norm : m_tangent_norm).block(0, 0, N, N);
    auto v_norm = m_v_norm.block(0, 0, N, N);
    auto h1 = m_h1.block(0, 0, N, N), h0 = m_h0.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto temp = h0 - h1 * ksqrtc * v_norm;
    auto vdotn = (transp? -1.0 : 1.0) * (v.real() * tangent.imag() - v.imag() * tangent.real());
    auto zmat = m_zero.block(0, 0, N, N);
    auto otherwise = (v_norm > epsilon).select(M_2_PI / (v_norm * v_norm), zmat);
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            auto fg = (transp ? m_double_layer_general_fg_t : m_double_layer_general_fg).block(i * N, j * N, N, N);
            auto cf = m_wg * t_norm * vdotn * fg;
            integral += (cf * mask.select(ksqrtc * h1 / v_norm, otherwise)).sum();
            if (der > 0)
                integral_der += ksqrtc * mask.select(h0 * cf, zmat).sum();
            if (der > 1)
                integral_der2 += mask.select(temp * cf, zmat).sum();
            double_layer_interaction_matrix(i, j) = integral;
            if (der > 0)
                double_layer_der_interaction_matrix(i, j) = sqrtc * integral_der;
            if (der > 1)
                double_layer_der2_interaction_matrix(i, j) = c * integral_der2;
        }
    }
}

// compute interaction matrices for the W submatrix, coinciding panels case
inline void GalerkinMatrixBuilder::hypersingular_coinciding(int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    auto tangent_p = m_tangent_p.block(0, 0, N, N), tangent = m_tangent.block(0, 0, N, N);
    auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
    auto h0 = m_h0.block(0, 0, N, N), h1 = m_h1.block(0, 0, N, N);
    auto v_norm = m_v_norm.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto otherwise = (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero);
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            if (i < j) { // symmetry
                hypersingular_interaction_matrix(i, j) = hypersingular_interaction_matrix(j, i);
                if (der > 0)
                    hypersingular_der_interaction_matrix(i, j) = hypersingular_der_interaction_matrix(j, i);
                if (der > 1)
                    hypersingular_der2_interaction_matrix(i, j) = hypersingular_der2_interaction_matrix(j, i);
                continue;
            }
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (K = 0; K < 2; ++K) {
                auto fg = m_hypersingular_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N);
                auto fg_arc = m_hypersingular_coinciding_fg_arc.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N);
                m_cf = fg_arc - kkc * (m_temp = fg * tdottp);
                integral += (m_wc * m_cf * mask.select(h0, otherwise)).sum();
                if (der > 0)
                    integral_der += mask.select(m_wc * (h1 * v_norm * m_cf + 2.0 * h0 * ksqrtc * m_temp), m_zero).sum();
                if (der > 1)
                    integral_der2 += mask.select(
                        m_wc * ((h0 * v_norm - h1 / ksqrtc) * m_cf * v_norm - 2.0 * m_temp * (2.0 * h1 * ksqrtc * v_norm - h0)), m_zero).sum();
            }
            hypersingular_interaction_matrix(i, j) = integral;
            if (der > 0)
                hypersingular_der_interaction_matrix(i, j) = -sqrtc * integral_der;
            if (der > 1)
                hypersingular_der2_interaction_matrix(i, j) = -c * integral_der2;
        }
    }
}

// compute interaction matrices for the W submatrix, adjacent panels case
inline void GalerkinMatrixBuilder::hypersingular_adjacent(bool swap, int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (K = 0; K < 2; ++K) {
                auto fg = swap ? (K > 0 ? m_hypersingular_adjacent_fg.block(j * N, i * N, N, N)
                                        : m_hypersingular_adjacent_fg_swap.block(i * N, j * N, N, N))
                               : (K > 0 ? m_hypersingular_adjacent_fg_swap.block(j * N, i * N, N, N)
                                        : m_hypersingular_adjacent_fg.block(i * N, j * N, N, N));
                auto fg_arc = swap ? (K > 0 ? m_hypersingular_adjacent_fg_arc.block(j * N, i * N, N, N)
                                            : m_hypersingular_adjacent_fg_arc_swap.block(i * N, j * N, N, N))
                                   : (K > 0 ? m_hypersingular_adjacent_fg_arc_swap.block(j * N, i * N, N, N)
                                            : m_hypersingular_adjacent_fg_arc.block(i * N, j * N, N, N));
                auto tangent_p = m_tangent_p.block(0, K * N, N, N), tangent = m_tangent.block(0, K * N, N, N);
                auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
                auto h0 = m_h0.block(0, K * N, N, N), h1 = m_h1.block(0, K * N, N, N);
                auto v_norm = m_v_norm.block(0, K * N, N, N);
                auto mask = ksqrtca * v_norm > epsilon;
                auto otherwise = (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero);
                m_cf = fg_arc - kkc * (m_temp = fg * tdottp);
                integral += (m_wa * m_cf * mask.select(h0, otherwise)).sum();
                if (der > 0)
                    integral_der += mask.select(m_wa * (h1 * v_norm * m_cf + 2.0 * h0 * ksqrtc * m_temp), m_zero).sum();
                if (der > 1)
                    integral_der2 += mask.select(
                        m_wa * ((h0 * v_norm - h1 / ksqrtc) * m_cf * v_norm - 2.0 * m_temp * (2.0 * h1 * ksqrtc * v_norm - h0)), m_zero).sum();
            }
            hypersingular_interaction_matrix(i, j) = integral;
            if (der > 0)
                hypersingular_der_interaction_matrix(i, j) = -sqrtc * integral_der;
            if (der > 1)
                hypersingular_der2_interaction_matrix(i, j) = -c * integral_der2;
        }
    }
}

// compute interaction matrices for the W submatrix, disjoint panels case
inline void GalerkinMatrixBuilder::hypersingular_general(int der) throw() {
    size_t i, j, N = GaussQR.n;
    auto tangent_p = m_tangent_p.block(0, 0, N, N), tangent = m_tangent.block(0, 0, N, N);
    auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
    auto h0 = m_h0.block(0, 0, N, N), h1 = m_h1.block(0, 0, N, N);
    auto v_norm = m_v_norm.block(0, 0, N, N);
    auto zmat = m_zero.block(0, 0, N, N);
    auto cf = m_cf.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto otherwise = (v_norm > epsilon).select(-M_2_PI * v_norm.log(), zmat);
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            auto fg = m_hypersingular_general_fg.block(i * N, j * N, N, N);
            auto fg_arc = m_hypersingular_general_fg_arc.block(i * N, j * N, N, N);
            cf = fg_arc - kkc * (m_temp = fg * tdottp);
            integral += (m_wg * cf * mask.select(h0, otherwise)).sum();
            if (der > 0)
                integral_der += mask.select(m_wg * (h1 * v_norm * cf + 2.0 * h0 * ksqrtc * m_temp), zmat).sum();
            if (der > 1)
                integral_der2 += mask.select(
                    m_wg * ((h0 * v_norm - h1 / ksqrtc) * cf * v_norm - 2.0 * m_temp * (2.0 * h1 * ksqrtc * v_norm - h0)), zmat).sum();
            hypersingular_interaction_matrix(i, j) = integral;
            if (der > 0)
                hypersingular_der_interaction_matrix(i, j) = -sqrtc * integral_der;
            if (der > 1)
                hypersingular_der2_interaction_matrix(i, j) = -c * integral_der2;
        }
    }
}

// compute interaction matrices for the V submatrix, coinciding panels case
inline void GalerkinMatrixBuilder::single_layer_coinciding(int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    auto t_norm = m_tangent_norm.block(0, 0, N, N);
    auto tp_norm = m_tangent_p_norm.block(0, 0, N, N);
    auto v_norm = m_v_norm.block(0, 0, N, N);
    auto h0 = m_h0.block(0, 0, N, N), h1 = m_h1.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto otherwise = (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero);
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            if (i < j) { // symmetry
                single_layer_interaction_matrix(i, j) = single_layer_interaction_matrix(j, i);
                if (der > 0)
                    single_layer_der_interaction_matrix(i, j) = single_layer_der_interaction_matrix(j, i);
                if (der > 1)
                    single_layer_der2_interaction_matrix(i, j) = single_layer_der2_interaction_matrix(j, i);
                continue;
            }
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (K = 0; K < 2; ++K) {
                auto fg = m_single_layer_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N);
                m_cf = tp_norm * t_norm * fg * m_wc;
                integral += (m_cf * mask.select(h0, otherwise)).sum();
                if (der > 0)
                    integral_der += mask.select(h1 * v_norm * m_cf, m_zero).sum();
                if (der > 1)
                    integral_der2 += mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm * m_cf, m_zero).sum();
            }
            single_layer_interaction_matrix(i, j) = integral;
            if (der > 0)
                single_layer_der_interaction_matrix(i, j) = -sqrtc * integral_der;
            if (der > 1)
                single_layer_der2_interaction_matrix(i, j) = c * integral_der2;
        }
    }
}

// compute interaction matrices for the V submatrix, adjacent panels case
inline void GalerkinMatrixBuilder::single_layer_adjacent(bool swap, int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (K = 0; K < 2; ++K) {
                auto fg = swap ? (K > 0 ? m_single_layer_adjacent_fg.block(j * N, i * N, N, N)
                                        : m_single_layer_adjacent_fg_swap.block(i * N, j * N, N, N))
                               : (K > 0 ? m_single_layer_adjacent_fg_swap.block(j * N, i * N, N, N)
                                        : m_single_layer_adjacent_fg.block(i * N, j * N, N, N));
                auto t_norm = m_tangent_norm.block(0, K * N, N, N);
                auto tp_norm = m_tangent_p_norm.block(0, K * N, N, N);
                auto v_norm = m_v_norm.block(0, K * N, N, N);
                auto h0 = m_h0.block(0, K * N, N, N), h1 = m_h1.block(0, K * N, N, N);
                auto mask = ksqrtca * v_norm > epsilon;
                auto otherwise = (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero);
                m_cf = tp_norm * t_norm * fg * m_wa;
                integral += (m_cf * mask.select(h0, otherwise)).sum();
                if (der > 0)
                    integral_der += mask.select(h1 * v_norm * m_cf, m_zero).sum();
                if (der > 1)
                    integral_der2 += mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm * m_cf, m_zero).sum();
            }
            single_layer_interaction_matrix(i, j) = integral;
            if (der > 0)
                single_layer_der_interaction_matrix(i, j) = -sqrtc * integral_der;
            if (der > 1)
                single_layer_der2_interaction_matrix(i, j) = c * integral_der2;
        }
    }
}

// compute interaction matrices for the V submatrix, disjoint panels case
inline void GalerkinMatrixBuilder::single_layer_general(int der) throw() {
    size_t i, j, N = GaussQR.n;
    auto t_norm = m_tangent_norm.block(0, 0, N, N);
    auto tp_norm = m_tangent_p_norm.block(0, 0, N, N);
    auto v_norm = m_v_norm.block(0, 0, N, N);
    auto h0 = m_h0.block(0, 0, N, N), h1 = m_h1.block(0, 0, N, N);
    auto zmat = m_zero.block(0, 0, N, N);
    auto cf = m_cf.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto otherwise = (v_norm > epsilon).select(-M_2_PI * v_norm.log(), zmat);
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            auto fg = m_single_layer_general_fg.block(i * N, j * N, N, N);
            cf = tp_norm * t_norm * fg * m_wg;
            integral += (cf * mask.select(h0, otherwise)).sum();
            if (der > 0)
                integral_der += mask.select(h1 * v_norm * cf, zmat).sum();
            if (der > 1)
                integral_der2 += mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm * cf, zmat).sum();
            single_layer_interaction_matrix(i, j) = integral;
            if (der > 0)
                single_layer_der_interaction_matrix(i, j) = -sqrtc * integral_der;
            if (der > 1)
                single_layer_der2_interaction_matrix(i, j) = c * integral_der2;
        }
    }
}

// initialize wavenumber and refraction index with related values
void GalerkinMatrixBuilder::initialize_parameters(const std::complex<double>& k_in, double c_in) {
    if (c_in < 1.)
        throw std::runtime_error("Refraction index must not be smaller than 1!");
    k = k_in;
    c = c_in;
    sqrtc = sqrt(c);
    ksqrtc = k * sqrtc;
    ksqrtca = std::abs(ksqrtc);
    kkc = ksqrtc * ksqrtc;
    k_real_positive = k.imag() == 0. && k.real() > 0;
}

// return true iff the panels P1 and P2 are adjacent, compute the SWAP boolean alongside
bool GalerkinMatrixBuilder::is_adjacent(const AbstractParametrizedCurve& p1, const AbstractParametrizedCurve& p2, bool &swap) const {
    double t1 = (p1(1) - p2(-1)).norm() / 100., t2 = (p1(-1) - p2(1)).norm() / 100.;
    if (t1 < epsilon || t2 < epsilon) {
        swap = t1 > epsilon;
        return true;
    }
    return false;
}

void GalerkinMatrixBuilder::assembleDoubleLayer(const std::complex<double>& k_in, double c_in, int der) {
    initialize_parameters(k_in, c_in);
    double_layer_matrix.setZero();
    if (der > 0)
        double_layer_der_matrix.setZero();
    if (der > 1)
        double_layer_der2_matrix.setZero();
    size_t i, j, I, J, II, JJ;
    bool swap;
    // Panel oriented assembly
    for (i = 0; i < dim_test; ++i) {
        const auto &pi = *panels[i];
        for (j = 0; j < dim_trial; ++j) {
            const auto &pj = *panels[j];
            if (i == j) {
                // coinciding panels
                compute_coinciding(pi);
                double_layer_coinciding(der);
            } else if (is_adjacent(pi, pj, swap)) {
                // adjacent panels
                compute_adjacent(pi, pj, swap);
                double_layer_adjacent(swap, der, false);
            } else {
                // disjoint panels
                compute_general(pi, pj);
                double_layer_general(der, false);
            }
            // Local to global mapping of the elements in interaction matrix
            for (I = 0; I < Qtest; ++I) {
                for (J = 0; J < Qtrial; ++J) {
                    II = test_space.LocGlobMap(I + 1, i + 1, dim_test) - 1;
                    JJ = trial_space.LocGlobMap(J + 1, j + 1, dim_trial) - 1;
                    double_layer_matrix(II, JJ) += double_layer_interaction_matrix(I, J);
                    if (der > 0)
                        double_layer_der_matrix(II, JJ) += double_layer_der_interaction_matrix(I, J);
                    if (der > 1)
                        double_layer_der2_matrix(II, JJ) += double_layer_der2_interaction_matrix(I, J);
                }
            }
        }
    }
    double_layer_matrix *= 0.25;
    if (der > 1)
        double_layer_der_matrix *= 0.25;
    if (der > 2)
        double_layer_der2_matrix *= 0.25;
}

void GalerkinMatrixBuilder::assembleHypersingular(const std::complex<double>& k_in, double c_in, int der) {
    initialize_parameters(k_in, c_in);
    hypersingular_matrix.setZero();
    if (der > 0)
        hypersingular_der_matrix.setZero();
    if (der > 1)
        hypersingular_der2_matrix.setZero();
    size_t i, j, I, J, II, JJ;
    bool swap;
    // Panel oriented assembly
    for (i = 0; i < dim_trial; ++i) {
        const auto &pi = *panels[i];
        for (j = 0; j <= i; ++j) {
            const auto &pj = *panels[j];
            if (i == j) {
                // coinciding panels
                compute_coinciding(pi);
                hypersingular_coinciding(der);
            } else if (is_adjacent(pi, pj, swap)) {
                // adjacent panels
                compute_adjacent(pi, pj, swap);
                hypersingular_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(pi, pj);
                hypersingular_general(der);
            }
            // Local to global mapping of the elements in interaction matrix
            for (I = 0; I < Qtrial; ++I) {
                for (J = 0; J < Qtrial; ++J) {
                    II = trial_space.LocGlobMap(I + 1, i + 1, dim_trial) - 1;
                    JJ = trial_space.LocGlobMap(J + 1, j + 1, dim_trial) - 1;
                    hypersingular_matrix(II, JJ) += hypersingular_interaction_matrix(I, J);
                    if (der > 0)
                        hypersingular_der_matrix(II, JJ) += hypersingular_der_interaction_matrix(I, J);
                    if (der > 1)
                        hypersingular_der2_matrix(II, JJ) += hypersingular_der2_interaction_matrix(I, J);
                }
            }
            if (i == j)
                continue;
            for (I = 0; I < Qtrial; ++I) {
                for (J = 0; J < Qtrial; ++J) {
                    II = trial_space.LocGlobMap(I + 1, j + 1, dim_trial) - 1;
                    JJ = trial_space.LocGlobMap(J + 1, i + 1, dim_trial) - 1;
                    hypersingular_matrix(II, JJ) += hypersingular_interaction_matrix(J, I);
                    if (der > 0)
                        hypersingular_der_matrix(II, JJ) += hypersingular_der_interaction_matrix(J, I);
                    if (der > 1)
                        hypersingular_der2_matrix(II, JJ) += hypersingular_der2_interaction_matrix(J, I);
                }
            }
        }
    }
    hypersingular_matrix *= 0.25;
    if (der > 1)
        hypersingular_der_matrix *= 0.25;
    if (der > 2)
        hypersingular_der2_matrix *= 0.25;
}

void GalerkinMatrixBuilder::assembleSingleLayer(const std::complex<double>& k_in, double c_in, int der) {
    initialize_parameters(k_in, c_in);
    single_layer_matrix.setZero();
    if (der > 0)
        single_layer_der_matrix.setZero();
    if (der > 1)
        single_layer_der2_matrix.setZero();
    size_t i, j, I, J, II, JJ;
    bool swap;
    // Panel oriented assembly
    for (i = 0; i < dim_test; ++i) {
        const auto &pi = *panels[i];
        for (j = 0; j <= i; ++j) {
            const auto &pj = *panels[j];
            if (i == j) {
                // coinciding panels
                compute_coinciding(pi);
                single_layer_coinciding(der);
            } else if (is_adjacent(pi, pj, swap)) {
                // adjacent panels
                compute_adjacent(pi, pj, swap);
                single_layer_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(pi, pj);
                single_layer_general(der);
            }
            // Local to global mapping of the elements in interaction matrix
            for (I = 0; I < Qtest; ++I) {
                for (J = 0; J < Qtest; ++J) {
                    II = test_space.LocGlobMap(I + 1, i + 1, dim_test) - 1;
                    JJ = test_space.LocGlobMap(J + 1, j + 1, dim_test) - 1;
                    single_layer_matrix(II, JJ) += single_layer_interaction_matrix(I, J);
                    if (der > 0)
                        single_layer_der_matrix(II, JJ) += single_layer_der_interaction_matrix(I, J);
                    if (der > 1)
                        single_layer_der2_matrix(II, JJ) += single_layer_der2_interaction_matrix(I, J);
                }
            }
            if (i == j)
                continue;
            for (I = 0; I < Qtest; ++I) {
                for (J = 0; J < Qtest; ++J) {
                    II = test_space.LocGlobMap(I + 1, j + 1, dim_test) - 1;
                    JJ = test_space.LocGlobMap(J + 1, i + 1, dim_test) - 1;
                    single_layer_matrix(II, JJ) += single_layer_interaction_matrix(J, I);
                    if (der > 0)
                        single_layer_der_matrix(II, JJ) += single_layer_der_interaction_matrix(J, I);
                    if (der > 1)
                        single_layer_der2_matrix(II, JJ) += single_layer_der2_interaction_matrix(J, I);
                }
            }
        }
    }
    single_layer_matrix *= 0.25;
    if (der > 1)
        single_layer_der_matrix *= 0.25;
    if (der > 2)
        single_layer_der2_matrix *= 0.25;
}

void GalerkinMatrixBuilder::assembleAll(const std::complex<double>& k_in, double c_in, int der) {
    if (&test_space != &trial_space)
        throw std::runtime_error("Trial and test spaces must be equal!");
    initialize_parameters(k_in, c_in);
    single_layer_matrix.setZero();
    double_layer_matrix.setZero();
    hypersingular_matrix.setZero();
    if (der > 0) {
        single_layer_der_matrix.setZero();
        double_layer_der_matrix.setZero();
        hypersingular_der_matrix.setZero();
    }
    if (der > 1) {
        single_layer_der2_matrix.setZero();
        double_layer_der2_matrix.setZero();
        hypersingular_der2_matrix.setZero();
    }
    size_t i, j, I, J, II, JJ, Q = Qtest, dim = dim_test;
    bool swap, adj;
    // Panel oriented assembly
    for (i = 0; i < dim; ++i) {
        const auto &pi = *panels[i];
        for (j = 0; j <= i; ++j) {
            const auto &pj = *panels[j];
            if (i == j) {
                // coinciding panels
                compute_coinciding(pi);
                auto tic = chrono::high_resolution_clock::now();
                double_layer_coinciding(der);
                hypersingular_coinciding(der);
                single_layer_coinciding(der);
                auto toc = chrono::high_resolution_clock::now();
                auto dur = chrono::duration_cast<chrono::microseconds>(toc - tic);
                interaction_matrix_time += dur.count();
            } else if ((adj = is_adjacent(pi, pj, swap))) {
                // adjacent panels
                compute_adjacent(pi, pj, swap);
                auto tic = chrono::high_resolution_clock::now();
                double_layer_adjacent(swap, der, false);
                hypersingular_adjacent(swap, der);
                single_layer_adjacent(swap, der);
                auto toc = chrono::high_resolution_clock::now();
                auto dur = chrono::duration_cast<chrono::microseconds>(toc - tic);
                interaction_matrix_time += dur.count();
            } else {
                // disjoint panels
                compute_general(pi, pj);
                auto tic = chrono::high_resolution_clock::now();
                double_layer_general(der, false);
                hypersingular_general(der);
                single_layer_general(der);
                auto toc = chrono::high_resolution_clock::now();
                auto dur = chrono::duration_cast<chrono::microseconds>(toc - tic);
                interaction_matrix_time += dur.count();
            }
            // Local to global mapping of the elements in interaction matrix
            auto tic = chrono::high_resolution_clock::now();
            for (I = 0; I < Q; ++I) {
                for (J = 0; J < Q; ++J) {
                    II = test_space.LocGlobMap(I + 1, i + 1, dim) - 1;
                    JJ = test_space.LocGlobMap(J + 1, j + 1, dim) - 1;
                    double_layer_matrix(II, JJ) += double_layer_interaction_matrix(I, J);
                    hypersingular_matrix(II, JJ) += hypersingular_interaction_matrix(I, J);
                    single_layer_matrix(II, JJ) += single_layer_interaction_matrix(I, J);
                    if (der > 0) {
                        double_layer_der_matrix(II, JJ) += double_layer_der_interaction_matrix(I, J);
                        hypersingular_der_matrix(II, JJ) += hypersingular_der_interaction_matrix(I, J);
                        single_layer_der_matrix(II, JJ) += single_layer_der_interaction_matrix(I, J);
                    }
                    if (der > 1) {
                        double_layer_der2_matrix(II, JJ) += double_layer_der2_interaction_matrix(I, J);
                        hypersingular_der2_matrix(II, JJ) += hypersingular_der2_interaction_matrix(I, J);
                        single_layer_der2_matrix(II, JJ) += single_layer_der2_interaction_matrix(I, J);
                    }
                }
            }
            auto toc = chrono::high_resolution_clock::now();
            auto dur = chrono::duration_cast<chrono::microseconds>(toc - tic);
            mapping_time += dur.count();
            if (i == j)
                continue;
            // use already computed data for the (j, i) case
            tic = chrono::high_resolution_clock::now();
            if (adj)
                double_layer_adjacent(!swap, der, true);
            else
                double_layer_general(der, true);
            toc = chrono::high_resolution_clock::now();
            dur = chrono::duration_cast<chrono::microseconds>(toc - tic);
            interaction_matrix_time += dur.count();
            tic = chrono::high_resolution_clock::now();
            for (I = 0; I < Q; ++I) {
                for (J = 0; J < Q; ++J) {
                    II = test_space.LocGlobMap(I + 1, j + 1, dim) - 1;
                    JJ = test_space.LocGlobMap(J + 1, i + 1, dim) - 1;
                    double_layer_matrix(II, JJ) += double_layer_interaction_matrix(I, J);
                    hypersingular_matrix(II, JJ) += hypersingular_interaction_matrix(J, I);
                    single_layer_matrix(II, JJ) += single_layer_interaction_matrix(J, I);
                    if (der > 0) {
                        double_layer_der_matrix(II, JJ) += double_layer_der_interaction_matrix(I, J);
                        hypersingular_der_matrix(II, JJ) += hypersingular_der_interaction_matrix(J, I);
                        single_layer_der_matrix(II, JJ) += single_layer_der_interaction_matrix(J, I);
                    }
                    if (der > 1) {
                        double_layer_der2_matrix(II, JJ) += double_layer_der2_interaction_matrix(I, J);
                        hypersingular_der2_matrix(II, JJ) += hypersingular_der2_interaction_matrix(J, I);
                        single_layer_der2_matrix(II, JJ) += single_layer_der2_interaction_matrix(J, I);
                    }
                }
            }
            toc = chrono::high_resolution_clock::now();
            dur = chrono::duration_cast<chrono::microseconds>(toc - tic);
            mapping_time += dur.count();
        }
    }
    double_layer_matrix *= 0.25;
    hypersingular_matrix *= 0.25;
    single_layer_matrix *= 0.25;
    if (der > 0) {
        double_layer_der_matrix *= 0.25;
        hypersingular_der_matrix *= 0.25;
        single_layer_der_matrix *= 0.25;
    }
    if (der > 1) {
        double_layer_der2_matrix *= 0.25;
        hypersingular_der2_matrix *= 0.25;
        single_layer_der2_matrix *= 0.25;
    }
}

const Eigen::MatrixXcd & GalerkinMatrixBuilder::getDoubleLayer(int der) const {
    if (der == 0)
        return double_layer_matrix;
    if (der == 1)
        return double_layer_der_matrix;
    if (der == 2)
        return double_layer_der2_matrix;
    throw std::runtime_error("Invalid order of derivative!");
}

const Eigen::MatrixXcd & GalerkinMatrixBuilder::getHypersingular(int der) const {
    if (der == 0)
        return hypersingular_matrix;
    if (der == 1)
        return hypersingular_der_matrix;
    if (der == 2)
        return hypersingular_der2_matrix;
    throw std::runtime_error("Invalid order of derivative!");
}

const Eigen::MatrixXcd & GalerkinMatrixBuilder::getSingleLayer(int der) const {
    if (der == 0)
        return single_layer_matrix;
    if (der == 1)
        return single_layer_der_matrix;
    if (der == 2)
        return single_layer_der2_matrix;
    throw std::runtime_error("Invalid order of derivative!");
}
