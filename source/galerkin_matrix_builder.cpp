#include "galerkin_matrix_builder.hpp"
#include "cbessel.hpp"
#include <numeric>
#include <execution>
#include <iostream>

using namespace std::complex_literals;

typedef std::complex<double> complex_t;
static const complex_t czero(0., 0.);
static const double epsilon = std::numeric_limits<double>::epsilon();

GalerkinMatrixBuilder::GalerkinMatrixBuilder(const ParametrizedMesh &mesh,
                                             const AbstractBEMSpace &test_space_in,
                                             const AbstractBEMSpace &trial_space_in,
                                             unsigned order)
: test_space(test_space_in), trial_space(trial_space_in)
{
    GaussQR = getGaussQR(order, 0., 1.);
    CGaussQR = getCGaussQR(order);
    panels = mesh.getPanels();
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
inline void GalerkinMatrixBuilder::compute_coinciding(const AbstractParametrizedCurve &p) throw() {
    size_t I, J, N = CGaussQR.n;
    m_tangent_p.block(0, 0, N, N) = p.Derivative_01(m_tc);
    m_tangent.block(0, 0, N, N) = p.Derivative_01(m_sc);
    m_v.block(0, 0, N, N) = p[m_sc] - p[m_tc];
    m_v_norm.block(0, 0, N, N) = m_v.block(0, 0, N, N).cwiseAbs();
    m_tangent_norm.block(0, 0, N, N) = m_tangent.block(0, 0, N, N).cwiseAbs();
    m_tangent_p_norm.block(0, 0, N, N) = m_tangent_p.block(0, 0, N, N).cwiseAbs();
    if (k_real_positive) {
        m_h_arg_half = ksqrtc.real() * m_v_norm.block(0, 0, N, N);
        complex_bessel::H1_01_i(m_h_arg_half, m_h0_res_half, m_h1_res_half);
        m_h0.block(0, 0, N, N) = m_h0_res_half;
        m_h1.block(0, 0, N, N) = m_h1_res_half;
    } else {
        for (J = 0; J < N; ++J) {
            for (I = 0; I < N; ++I) {
                m_h0(I, J) = 1i * complex_bessel::H1(0, ksqrtc * m_v_norm(I, J));
                m_h1(I, J) = 1i * complex_bessel::H1(1, ksqrtc * m_v_norm(I, J));
            }
        }
    }
}

// compute values required for the adjacent panels case
inline void GalerkinMatrixBuilder::compute_adjacent(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p, bool swap) throw() {
    size_t I, J, N = CGaussQR.n;
    m_tangent_p.block(0, 0, N, N) = swap ? -pi_p.Derivative_01_swapped(m_ta) : pi_p.Derivative_01(m_ta);
    m_tangent_p.block(0, N, N, N) = swap ? -pi_p.Derivative_01_swapped(m_sa) : pi_p.Derivative_01(m_sa);
    m_tangent.block(0, 0, N, N) = swap ? pi.Derivative_01(m_sa) : -pi.Derivative_01_swapped(m_sa);
    m_tangent.block(0, N, N, N) = swap ? pi.Derivative_01(m_ta) : -pi.Derivative_01_swapped(m_ta);
    m_v.block(0, 0, N, N) = swap ? pi[m_sa] - pi_p.swapped_op(m_ta) : pi.swapped_op(m_sa) - pi_p[m_ta];
    m_v.block(0, N, N, N) = swap ? pi[m_ta] - pi_p.swapped_op(m_sa) : pi.swapped_op(m_ta) - pi_p[m_sa];
    m_v_norm = m_v.cwiseAbs();
    m_tangent_norm = m_tangent.cwiseAbs();
    m_tangent_p_norm = m_tangent_p.cwiseAbs();
    if (k_real_positive) {
        m_h_arg = ksqrtc.real() * m_v_norm;
        complex_bessel::H1_01_i(m_h_arg, m_h0, m_h1);
    } else {
        for (J = 0; J < 2 * N; ++J) {
            for (I = 0; I < N; ++I) {
                m_h0(I, J) = 1i * complex_bessel::H1(0, ksqrtc * m_v_norm(I, J));
                m_h1(I, J) = 1i * complex_bessel::H1(1, ksqrtc * m_v_norm(I, J));
            }
        }
    }
}

// compute values required for the disjoint panels case
inline void GalerkinMatrixBuilder::compute_general(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p) throw() {
    size_t I, J, N = GaussQR.n;
    m_tangent_p.block(0, 0, N, N) = pi_p.Derivative_01(m_tg);
    m_tangent.block(0, 0, N, N) = pi.Derivative_01(m_sg);
    m_v.block(0, 0, N, N) = pi[m_sg] - pi_p[m_tg];
    m_v_norm.block(0, 0, N, N) = m_v.block(0, 0, N, N).cwiseAbs();
    m_tangent_norm.block(0, 0, N, N) = m_tangent.block(0, 0, N, N).cwiseAbs();
    m_tangent_p_norm.block(0, 0, N, N) = m_tangent_p.block(0, 0, N, N).cwiseAbs();
    if (k_real_positive) {
        m_h_arg_small = ksqrtc.real() * m_v_norm.block(0, 0, N, N);
        complex_bessel::H1_01_i(m_h_arg_small, m_h0_res_small, m_h1_res_small);
        m_h0.block(0, 0, N, N) = m_h0_res_small;
        m_h1.block(0, 0, N, N) = m_h1_res_small;
    } else {
        for (J = 0; J < N; ++J) {
            for (I = 0; I < N; ++I) {
                m_h0(I, J) = 1i * complex_bessel::H1(0, ksqrtc * m_v_norm(I, J));
                m_h1(I, J) = 1i * complex_bessel::H1(1, ksqrtc * m_v_norm(I, J));
            }
        }
    }
}

// compute interaction matrices for the K submatrix, coinciding panels case
void GalerkinMatrixBuilder::double_layer_coinciding(int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    const auto &v = m_v.block(0, 0, N, N);
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    const auto &h1 = m_h1.block(0, 0, N, N), &h0 = m_h0.block(0, 0, N, N);
    auto mask = (ksqrtca * v_norm) > epsilon;
    auto masked1 = mask.select(ksqrtc * h1 / v_norm, (v_norm > epsilon).select(M_2_PI / (v_norm * v_norm), m_zero));
    Eigen::ArrayXXcd masked2, masked3;
    double_layer_interaction_matrix.setZero();
    if (der > 0) {
        double_layer_der_interaction_matrix.setZero();
        masked2 = mask.select(h0, m_zero);
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix.setZero();
        masked3 = masked2 - mask.select(h1 * ksqrtc * v_norm, m_zero);
    }
    for (K = 0; K < 2; ++K) {
        const auto tangent = (K == 0 ? m_tangent_p : m_tangent).block(0, 0, N, N);
        m_temp = (K * 2 - 1.) * m_wc * (K == 0 ? m_tangent_norm : m_tangent_p_norm).block(0, 0, N, N) *
            (v.imag() * tangent.real() - v.real() * tangent.imag());
        for (j = 0; j < Qtrial; ++j) {
            for (i = 0; i < Qtest; ++i) {
                m_cf = m_temp * (K == 0 ? m_double_layer_coinciding_fg : m_double_layer_coinciding_fg_t).block(i * N, j * N, N, N);
                double_layer_interaction_matrix(i, j) += (m_cf * masked1).sum();
                if (der > 0)
                    double_layer_der_interaction_matrix(i, j) += (m_cf * masked2).sum();
                if (der > 1)
                    double_layer_der2_interaction_matrix(i, j) += (m_cf * masked3).sum();
            }
        }
    }
    if (der > 0)
        double_layer_der_interaction_matrix *= k * c;
    if (der > 1)
        double_layer_der2_interaction_matrix *= c;
}

// compute interaction matrices for the K submatrix, adjacent panels case
void GalerkinMatrixBuilder::double_layer_adjacent(bool swap, int der, bool transp) throw() {
    size_t i, j, K, N = CGaussQR.n;
    double_layer_interaction_matrix.setZero();
    if (der > 0)
        double_layer_der_interaction_matrix.setZero();
    if (der > 1)
        double_layer_der2_interaction_matrix.setZero();
    Eigen::ArrayXXcd masked2, masked3;
    for (K = 0; K < 2; ++K) {
        const auto &tangent = (transp ? m_tangent : m_tangent_p).block(0, K * N, N, N);
        const auto &v = m_v.block(0, K * N, N, N);
        const auto &v_norm = m_v_norm.block(0, K * N, N, N);
        const auto &h1 = m_h1.block(0, K * N, N, N), &h0 = m_h0.block(0, K * N, N, N);
        auto mask = (ksqrtca * v_norm) > epsilon;
        auto masked1 = mask.select(ksqrtc * h1 / v_norm, (v_norm > epsilon).select(M_2_PI / (v_norm * v_norm), m_zero));
        if (der > 0)
            masked2 = mask.select(h0, m_zero);
        if (der > 1)
            masked3 = masked2 - mask.select(h1 * ksqrtc * v_norm, m_zero);
        m_temp = (transp? 1.0 : -1.0) * m_wa * (v.imag() * tangent.real() - v.real() * tangent.imag()) *
            (transp ? m_tangent_p_norm : m_tangent_norm).block(0, K * N, N, N);
        for (j = 0; j < Qtrial; ++j) {
            for (i = 0; i < Qtest; ++i) {
                m_cf = m_temp * (swap ? ((transp ? K == 0 : K == 1) ? m_double_layer_adjacent_fg_swap_t
                                                                    : m_double_layer_adjacent_fg_swap)
                                      : ((transp ? K == 0 : K == 1) ? m_double_layer_adjacent_fg_t
                                                                    : m_double_layer_adjacent_fg)).block(i * N, j * N, N, N);
                double_layer_interaction_matrix(i, j) += (m_cf * masked1).sum();
                if (der > 0)
                    double_layer_der_interaction_matrix(i, j) += (m_cf * masked2).sum();
                if (der > 1)
                    double_layer_der2_interaction_matrix(i, j) += (m_cf * masked3).sum();
            }
        }
    }
    if (der > 0)
        double_layer_der_interaction_matrix *= k * c;
    if (der > 1)
        double_layer_der2_interaction_matrix *= c;
}

// compute interaction matrices for the K submatrix, disjoint panels case
void GalerkinMatrixBuilder::double_layer_general(int der, bool transp) throw() {
    size_t i, j, N = GaussQR.n;
    const auto &tangent = (transp ? m_tangent : m_tangent_p).block(0, 0, N, N), &v = m_v.block(0, 0, N, N);
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    const auto &h1 = m_h1.block(0, 0, N, N), &h0 = m_h0.block(0, 0, N, N);
    const auto &zmat = m_zero.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto masked1 = mask.select(ksqrtc * h1 / v_norm, (v_norm > epsilon).select(M_2_PI / (v_norm * v_norm), zmat));
    m_vdotn.block(0, 0, N, N) = m_wg * (transp? -1.0 : 1.0) * (v.real() * tangent.imag() - v.imag() * tangent.real()) *
        (transp ? m_tangent_p_norm : m_tangent_norm).block(0, 0, N, N);
    double_layer_interaction_matrix.setZero();
    Eigen::ArrayXXcd masked2, masked3;
    if (der > 0) {
        double_layer_der_interaction_matrix.setZero();
        masked2 = mask.select(h0, zmat);
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix.setZero();
        masked3 = masked2 - mask.select(h1 * ksqrtc * v_norm, zmat);
    }
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            m_cf.block(0, 0, N, N) = m_vdotn.block(0, 0, N, N) *
                (transp ? m_double_layer_general_fg_t : m_double_layer_general_fg).block(i * N, j * N, N, N);
            double_layer_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked1).sum();
            if (der > 0)
                double_layer_der_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked2).sum();
            if (der > 1)
                double_layer_der2_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked3).sum();
        }
    }
    if (der > 0)
        double_layer_der_interaction_matrix *= k * c;
    if (der > 1)
        double_layer_der2_interaction_matrix *= c;
}

// compute interaction matrices for the W submatrix, coinciding panels case
void GalerkinMatrixBuilder::hypersingular_coinciding(int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    const auto &tangent_p = m_tangent_p.block(0, 0, N, N), &tangent = m_tangent.block(0, 0, N, N);
    const auto &h0 = m_h0.block(0, 0, N, N), &h1 = m_h1.block(0, 0, N, N);
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
    auto mask = ksqrtca * v_norm > epsilon;
    auto masked1 = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
    hypersingular_interaction_matrix.setZero();
    if (der > 0)
        hypersingular_der_interaction_matrix.setZero();
    if (der > 1)
        hypersingular_der2_interaction_matrix.setZero();
    for (K = 0; K < 2; ++K) {
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
                m_cf = m_hypersingular_coinciding_fg_arc.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) -
                    kkc * (m_temp = m_hypersingular_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) * tdottp);
                hypersingular_interaction_matrix(i, j) += (m_wc * m_cf * masked1).sum();
                if (der > 0)
                    hypersingular_der_interaction_matrix(i, j) += mask.select(
                        m_wc * (h1 * v_norm * m_cf + 2.0 * h0 * ksqrtc * m_temp), m_zero).sum();
                if (der > 1)
                    hypersingular_der2_interaction_matrix(i, j) += mask.select(
                        m_wc * ((h0 * v_norm - h1 / ksqrtc) * m_cf * v_norm - 2.0 * m_temp * (2.0 * h1 * ksqrtc * v_norm - h0)), m_zero).sum();
            }
        }
    }
    if (der > 0)
        hypersingular_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        hypersingular_der2_interaction_matrix *= -c;
}

// compute interaction matrices for the W submatrix, adjacent panels case
void GalerkinMatrixBuilder::hypersingular_adjacent(bool swap, int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    hypersingular_interaction_matrix.setZero();
    if (der > 0)
        hypersingular_der_interaction_matrix.setZero();
    if (der > 1)
        hypersingular_der2_interaction_matrix.setZero();
    for (K = 0; K < 2; ++K) {
        const auto &tangent_p = m_tangent_p.block(0, K * N, N, N), &tangent = m_tangent.block(0, K * N, N, N);
        const auto &h0 = m_h0.block(0, K * N, N, N), &h1 = m_h1.block(0, K * N, N, N);
        const auto &v_norm = m_v_norm.block(0, K * N, N, N);
        auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
        auto mask = ksqrtca * v_norm > epsilon;
        auto masked1 = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
        for (j = 0; j < Qtrial; ++j) {
            for (i = 0; i < Qtrial; ++i) {
                m_cf = (swap ? (K > 0 ? m_hypersingular_adjacent_fg_arc.block(j * N, i * N, N, N)
                                      : m_hypersingular_adjacent_fg_arc_swap.block(i * N, j * N, N, N))
                             : (K > 0 ? m_hypersingular_adjacent_fg_arc_swap.block(j * N, i * N, N, N)
                                      : m_hypersingular_adjacent_fg_arc.block(i * N, j * N, N, N))) -
                    kkc * (m_temp = (swap ? (K > 0 ? m_hypersingular_adjacent_fg.block(j * N, i * N, N, N)
                                                   : m_hypersingular_adjacent_fg_swap.block(i * N, j * N, N, N))
                                          : (K > 0 ? m_hypersingular_adjacent_fg_swap.block(j * N, i * N, N, N)
                                                   : m_hypersingular_adjacent_fg.block(i * N, j * N, N, N))) * tdottp);
                hypersingular_interaction_matrix(i, j) += (m_wa * m_cf * masked1).sum();
                if (der > 0)
                    hypersingular_der_interaction_matrix(i, j) += mask.select(
                        m_wa * (h1 * v_norm * m_cf + 2.0 * h0 * ksqrtc * m_temp), m_zero).sum();
                if (der > 1)
                    hypersingular_der2_interaction_matrix(i, j) += mask.select(
                        m_wa * ((h0 * v_norm - h1 / ksqrtc) * m_cf * v_norm - 2.0 * m_temp * (2.0 * h1 * ksqrtc * v_norm - h0)), m_zero).sum();
            }
        }
    }
    if (der > 0)
        hypersingular_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        hypersingular_der2_interaction_matrix *= -c;
}

// compute interaction matrices for the W submatrix, disjoint panels case
void GalerkinMatrixBuilder::hypersingular_general(int der) throw() {
    size_t i, j, N = GaussQR.n;
    const auto &tangent_p = m_tangent_p.block(0, 0, N, N), &tangent = m_tangent.block(0, 0, N, N);
    const auto &h0 = m_h0.block(0, 0, N, N), &h1 = m_h1.block(0, 0, N, N);
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    const auto &zmat = m_zero.block(0, 0, N, N);
    auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
    auto mask = ksqrtca * v_norm > epsilon;
    auto masked1 = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), zmat));
    hypersingular_interaction_matrix.setZero();
    if (der > 0)
        hypersingular_der_interaction_matrix.setZero();
    if (der > 1)
        hypersingular_der2_interaction_matrix.setZero();
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            m_cf.block(0, 0, N, N) = m_hypersingular_general_fg_arc.block(i * N, j * N, N, N) -
                kkc * (m_temp = m_hypersingular_general_fg.block(i * N, j * N, N, N) * tdottp);
            hypersingular_interaction_matrix(i, j) += (m_wg * m_cf.block(0, 0, N, N) * masked1).sum();
            if (der > 0)
                hypersingular_der_interaction_matrix(i, j) +=
                    (m_wg * mask.select(h1 * v_norm * m_cf.block(0, 0, N, N) + 2.0 * h0 * ksqrtc * m_temp, zmat)).sum();
            if (der > 1)
                hypersingular_der2_interaction_matrix(i, j) += (m_wg * mask.select(
                    (h0 * v_norm - h1 / ksqrtc) * m_cf.block(0, 0, N, N) * v_norm - 2.0 * m_temp * (2.0 * h1 * ksqrtc * v_norm - h0), zmat)).sum();
        }
    }
    if (der > 0)
        hypersingular_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        hypersingular_der2_interaction_matrix *= -c;
}

// compute interaction matrices for the V submatrix, coinciding panels case
void GalerkinMatrixBuilder::single_layer_coinciding(int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    const auto &h0 = m_h0.block(0, 0, N, N), &h1 = m_h1.block(0, 0, N, N);
    auto ttp_norm = m_tangent_norm.block(0, 0, N, N) * m_tangent_p_norm.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto masked1 = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
    Eigen::ArrayXXcd masked2, masked3;
    single_layer_interaction_matrix.setZero();
    if (der > 0) {
        single_layer_der_interaction_matrix.setZero();
        masked2 = mask.select(h1 * v_norm, m_zero);
    }
    if (der > 1) {
        single_layer_der2_interaction_matrix.setZero();
        masked3 = mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm, m_zero);
    }
    for (K = 0; K < 2; ++K) {
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
                m_cf = ttp_norm * m_single_layer_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) * m_wc;
                single_layer_interaction_matrix(i, j) += (m_cf * masked1).sum();
                if (der > 0)
                    single_layer_der_interaction_matrix(i, j) += (m_cf * masked2).sum();
                if (der > 1)
                    single_layer_der2_interaction_matrix(i, j) += (m_cf * masked3).sum();
            }
        }
    }
    if (der > 0)
        single_layer_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        single_layer_der2_interaction_matrix *= c;
}

// compute interaction matrices for the V submatrix, adjacent panels case
void GalerkinMatrixBuilder::single_layer_adjacent(bool swap, int der) throw() {
    size_t i, j, K, N = CGaussQR.n;
    auto ttp_norm_0 = m_tangent_norm.block(0, 0, N, N) * m_tangent_p_norm.block(0, 0, N, N);
    auto ttp_norm_1 = m_tangent_norm.block(0, N, N, N) * m_tangent_p_norm.block(0, N, N, N);
    single_layer_interaction_matrix.setZero();
    if (der > 0)
        single_layer_der_interaction_matrix.setZero();
    if (der > 1)
        single_layer_der2_interaction_matrix.setZero();
    Eigen::ArrayXXcd masked2, masked3;
    for (K = 0; K < 2; ++K) {
        const auto &v_norm = m_v_norm.block(0, K * N, N, N);
        const auto &h0 = m_h0.block(0, K * N, N, N), &h1 = m_h1.block(0, K * N, N, N);
        auto mask = ksqrtca * v_norm > epsilon;
        auto masked1 = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
        if (der > 0)
            masked2 = mask.select(h1 * v_norm, m_zero);
        if (der > 1)
            masked3 = mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm, m_zero);
        m_temp = (K == 0 ? ttp_norm_0 : ttp_norm_1) * m_wa;
        for (j = 0; j < Qtest; ++j) {
            for (i = 0; i < Qtest; ++i) {
                m_cf = m_temp * (swap ? (K > 0 ? m_single_layer_adjacent_fg.block(j * N, i * N, N, N)
                                               : m_single_layer_adjacent_fg_swap.block(i * N, j * N, N, N))
                                      : (K > 0 ? m_single_layer_adjacent_fg_swap.block(j * N, i * N, N, N)
                                               : m_single_layer_adjacent_fg.block(i * N, j * N, N, N)));
                single_layer_interaction_matrix(i, j) += (m_cf * masked1).sum();
                if (der > 0)
                    single_layer_der_interaction_matrix(i, j) += (m_cf * masked2).sum();
                if (der > 1)
                    single_layer_der2_interaction_matrix(i, j) += (m_cf * masked3).sum();
            }
        }
    }
    if (der > 0)
        single_layer_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        single_layer_der2_interaction_matrix *= c;
}

// compute interaction matrices for the V submatrix, disjoint panels case
void GalerkinMatrixBuilder::single_layer_general(int der) throw() {
    size_t i, j, N = GaussQR.n;
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    const auto &h0 = m_h0.block(0, 0, N, N), &h1 = m_h1.block(0, 0, N, N);
    const auto &zmat = m_zero.block(0, 0, N, N);
    auto ttp_norm = m_tangent_norm.block(0, 0, N, N) * m_tangent_p_norm.block(0, 0, N, N);
    auto mask = ksqrtca * v_norm > epsilon;
    auto masked1 = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), zmat));
    Eigen::ArrayXXcd masked2, masked3;
    single_layer_interaction_matrix.setZero();
    if (der > 0) {
        single_layer_der_interaction_matrix.setZero();
        masked2 = mask.select(h1 * v_norm, zmat);
    }
    if (der > 1) {
        single_layer_der2_interaction_matrix.setZero();
        masked3 = mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm, zmat);
    }
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            m_cf.block(0, 0, N, N) = ttp_norm * m_wg * m_single_layer_general_fg.block(i * N, j * N, N, N);
            single_layer_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked1).sum();
            if (der > 0)
                single_layer_der_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked2).sum();
            if (der > 1)
                single_layer_der2_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked3).sum();
        }
    }
    if (der > 0)
        single_layer_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        single_layer_der2_interaction_matrix *= c;
}

void GalerkinMatrixBuilder::all_coinciding(int der) throw() {
    size_t i, j, K, N = CGaussQR.n, Q = Qtest;
    const auto &v = m_v.block(0, 0, N, N);
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    const auto &h1 = m_h1.block(0, 0, N, N), &h0 = m_h0.block(0, 0, N, N);
    const auto &tangent_p = m_tangent_p.block(0, 0, N, N), &tangent = m_tangent.block(0, 0, N, N);
    auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
    auto ttp_norm = m_tangent_norm.block(0, 0, N, N) * m_tangent_p_norm.block(0, 0, N, N);
    auto v_norm2 = v_norm * v_norm;
    auto mask = (ksqrtca * v_norm) > epsilon;
    auto mask2 = v_norm > epsilon;
    auto masked1 = mask.select(ksqrtc * h1 / v_norm, mask2.select(M_2_PI / v_norm2, m_zero));
    auto masked1_hg = mask.select(h0, mask2.select(-M_2_PI * v_norm.log(), m_zero));
    Eigen::ArrayXXcd masked2, masked2_g, masked3, masked3_g;
    double_layer_interaction_matrix.setZero();
    hypersingular_interaction_matrix.setZero();
    single_layer_interaction_matrix.setZero();
    if (der > 0) {
        double_layer_der_interaction_matrix.setZero();
        hypersingular_der_interaction_matrix.setZero();
        single_layer_der_interaction_matrix.setZero();
        masked2 = mask.select(h0, m_zero);
        masked2_g = mask.select(h1 * v_norm, m_zero);
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix.setZero();
        hypersingular_der2_interaction_matrix.setZero();
        single_layer_der2_interaction_matrix.setZero();
        masked3 = masked2 - ksqrtc * masked2_g;
        masked3_g = masked2_g / ksqrtc - masked2 * v_norm2;
    }
    for (K = 0; K < 2; ++K) {
        auto temp = (K * 2 - 1.) * m_wc * (K == 0 ? m_tangent_norm : m_tangent_p_norm).block(0, 0, N, N) *
            (v.imag() * (K == 0 ? tangent_p : tangent).real() - v.real() * (K == 0 ? tangent_p : tangent).imag());
        for (j = 0; j < Q; ++j) {
            for (i = 0; i < Q; ++i) {
                m_cf = temp * (K == 0 ? m_double_layer_coinciding_fg : m_double_layer_coinciding_fg_t).block(i * N, j * N, N, N);
                double_layer_interaction_matrix(i, j) += (m_cf * masked1).sum();
                if (der > 0)
                    double_layer_der_interaction_matrix(i, j) += (m_cf * masked2).sum();
                if (der > 1)
                    double_layer_der2_interaction_matrix(i, j) += (m_cf * masked3).sum();
                if (i < j) { // symmetry
                    hypersingular_interaction_matrix(i, j) = hypersingular_interaction_matrix(j, i);
                    single_layer_interaction_matrix(i, j) = single_layer_interaction_matrix(j, i);
                    if (der > 0) {
                        hypersingular_der_interaction_matrix(i, j) = hypersingular_der_interaction_matrix(j, i);
                        single_layer_der_interaction_matrix(i, j) = single_layer_der_interaction_matrix(j, i);
                    }
                    if (der > 1) {
                        hypersingular_der2_interaction_matrix(i, j) = hypersingular_der2_interaction_matrix(j, i);
                        single_layer_der2_interaction_matrix(i, j) = single_layer_der2_interaction_matrix(j, i);
                    }
                } else {
                    m_cf = m_hypersingular_coinciding_fg_arc.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) -
                        kkc * (m_temp = m_hypersingular_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) * tdottp);
                    hypersingular_interaction_matrix(i, j) += (m_wc * m_cf * masked1_hg).sum();
                    if (der > 0)
                        hypersingular_der_interaction_matrix(i, j) += (m_wc * (masked2_g * m_cf + 2.0 * masked2 * ksqrtc * m_temp)).sum();
                    if (der > 1)
                        hypersingular_der2_interaction_matrix(i, j) +=
                            (m_wc * ((masked2 * v_norm2 - masked2_g / ksqrtc) * m_cf - 2.0 * m_temp * (masked2 - 2.0 * masked3))).sum();
                    m_cf = ttp_norm * m_single_layer_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) * m_wc;
                    single_layer_interaction_matrix(i, j) += (m_cf * masked1_hg).sum();
                    if (der > 0)
                        single_layer_der_interaction_matrix(i, j) += (m_cf * masked2_g).sum();
                    if (der > 1)
                        single_layer_der2_interaction_matrix(i, j) += (m_cf * masked3_g).sum();
                }
            }
        }
    }
    if (der > 0) {
        double_layer_der_interaction_matrix *= k * c;
        hypersingular_der_interaction_matrix *= -sqrtc;
        single_layer_der_interaction_matrix *= -sqrtc;
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix *= c;
        hypersingular_der2_interaction_matrix *= -c;
        single_layer_der2_interaction_matrix *= c;
    }
}

void GalerkinMatrixBuilder::all_general(int der) throw() {
    size_t i, j, N = GaussQR.n, Q = Qtest;
    const auto &tangent_p = m_tangent_p.block(0, 0, N, N), &tangent = m_tangent.block(0, 0, N, N);
    const auto &v = m_v.block(0, 0, N, N);
    const auto &v_norm = m_v_norm.block(0, 0, N, N);
    const auto &h1 = m_h1.block(0, 0, N, N), h0 = m_h0.block(0, 0, N, N);
    const auto &zmat = m_zero.block(0, 0, N, N);
    auto ttp_norm = m_tangent_norm.block(0, 0, N, N) * m_tangent_p_norm.block(0, 0, N, N);
    auto tdottp = tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
    auto v_norm2 = v_norm * v_norm;
    auto mask = ksqrtca * v_norm > epsilon;
    auto mask2 = v_norm > epsilon;
    auto masked1 = mask.select(ksqrtc * h1 / v_norm, mask2.select(M_2_PI / v_norm2, zmat));
    auto masked1_hg = mask.select(h0, mask2.select(-M_2_PI * v_norm.log(), zmat));
    m_vdotn.block(0, 0, N, N) = m_wg * (v.real() * tangent_p.imag() - v.imag() * tangent_p.real()) * m_tangent_norm.block(0, 0, N, N);
    Eigen::ArrayXXcd masked2, masked2_g, masked3, masked3_g;
    double_layer_interaction_matrix.setZero();
    hypersingular_interaction_matrix.setZero();
    single_layer_interaction_matrix.setZero();
    if (der > 0) {
        double_layer_der_interaction_matrix.setZero();
        hypersingular_der_interaction_matrix.setZero();
        single_layer_der_interaction_matrix.setZero();
        masked2 = mask.select(h0, zmat);
        masked2_g = mask.select(h1 * v_norm, zmat);
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix.setZero();
        hypersingular_der2_interaction_matrix.setZero();
        single_layer_der2_interaction_matrix.setZero();
        masked3 = masked2 - ksqrtc * masked2_g;
        masked3_g = masked2_g / ksqrtc - masked2 * v_norm2;
    }
    for (j = 0; j < Q; ++j) {
        for (i = 0; i < Q; ++i) {
            m_cf.block(0, 0, N, N) = m_vdotn.block(0, 0, N, N) * m_double_layer_general_fg.block(i * N, j * N, N, N);
            double_layer_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked1).sum();
            if (der > 0)
                double_layer_der_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked2).sum();
            if (der > 1)
                double_layer_der2_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked3).sum();
            m_cf.block(0, 0, N, N) = m_hypersingular_general_fg_arc.block(i * N, j * N, N, N) -
                kkc * (m_temp = m_hypersingular_general_fg.block(i * N, j * N, N, N) * tdottp);
            hypersingular_interaction_matrix(i, j) += (m_wg * m_cf.block(0, 0, N, N) * masked1_hg).sum();
            if (der > 0)
                hypersingular_der_interaction_matrix(i, j) +=
                    (m_wg * mask.select(h1 * v_norm * m_cf.block(0, 0, N, N) + 2.0 * h0 * ksqrtc * m_temp, zmat)).sum();
            if (der > 1)
                hypersingular_der2_interaction_matrix(i, j) += (m_wg * mask.select(
                    (h0 * v_norm - h1 / ksqrtc) * m_cf.block(0, 0, N, N) * v_norm - 2.0 * m_temp * (2.0 * h1 * ksqrtc * v_norm - h0), zmat)).sum();
            m_cf.block(0, 0, N, N) = ttp_norm * m_wg * m_single_layer_general_fg.block(i * N, j * N, N, N);
            single_layer_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked1_hg).sum();
            if (der > 0)
                single_layer_der_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked2_g).sum();
            if (der > 1)
                single_layer_der2_interaction_matrix(i, j) += (m_cf.block(0, 0, N, N) * masked3_g).sum();
        }
    }
    if (der > 0) {
        double_layer_der_interaction_matrix *= k * c;
        hypersingular_der_interaction_matrix *= -sqrtc;
        single_layer_der_interaction_matrix *= -sqrtc;
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix *= c;
        hypersingular_der2_interaction_matrix *= -c;
        single_layer_der2_interaction_matrix *= c;
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
            if (&pi == &pj) {
                // coinciding panels
                compute_coinciding(pi);
                all_coinciding(der);
            } else if ((adj = is_adjacent(pi, pj, swap))) {
                // adjacent panels
                compute_adjacent(pi, pj, swap);
                double_layer_adjacent(swap, der, false);
                hypersingular_adjacent(swap, der);
                single_layer_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(pi, pj);
                all_general(der);
            }
            // Local to global mapping of the elements in interaction matrix
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
            if (i == j)
                continue;
            // use already computed data for the (j, i) case
            if (adj)
                double_layer_adjacent(!swap, der, true);
            else
                double_layer_general(der, true);
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
