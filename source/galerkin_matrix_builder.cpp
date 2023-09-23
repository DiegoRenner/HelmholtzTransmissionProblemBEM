#include "galerkin_matrix_builder.hpp"
#include "cbessel.hpp"
#include <numeric>
#include <iostream>
#include <chrono>

static const double epsilon = std::numeric_limits<double>::epsilon();

BuilderData::BuilderData(const ParametrizedMesh &mesh_in,
                         const AbstractBEMSpace &test_space_in,
                         const AbstractBEMSpace &trial_space_in,
                         unsigned order)
: mesh(mesh_in), test_space(test_space_in), trial_space(trial_space_in)
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
    Derivative_01_sc.resize(N, N * numpanels);
    Derivative_01_tc.resize(N, N * numpanels);
    Derivative_01_sa.resize(N, N * numpanels);
    Derivative_01_ta.resize(N, N * numpanels);
    Derivative_01_swapped_sa.resize(N, N * numpanels);
    Derivative_01_swapped_ta.resize(N, N * numpanels);
    Derivative_01_sg.resize(Ns, Ns * numpanels);
    Derivative_01_tg.resize(Ns, Ns * numpanels);
    Derivative_01_sc_n.resize(N, N * numpanels);
    Derivative_01_tc_n.resize(N, N * numpanels);
    Derivative_01_sa_n.resize(N, N * numpanels);
    Derivative_01_ta_n.resize(N, N * numpanels);
    Derivative_01_swapped_sa_n.resize(N, N * numpanels);
    Derivative_01_swapped_ta_n.resize(N, N * numpanels);
    Derivative_01_sg_n.resize(Ns, Ns * numpanels);
    Derivative_01_tg_n.resize(Ns, Ns * numpanels);
    op_sc.resize(N, N * numpanels);
    op_tc.resize(N, N * numpanels);
    op_sa.resize(N, N * numpanels);
    op_ta.resize(N, N * numpanels);
    swapped_op_sa.resize(N, N * numpanels);
    swapped_op_ta.resize(N, N * numpanels);
    op_sg.resize(Ns, Ns * numpanels);
    op_tg.resize(Ns, Ns * numpanels);
    Eigen::ArrayXXcd tmp(N, N), tmps(Ns, Ns);
    Eigen::ArrayXXd tmp_n(N, N), tmps_n(Ns, Ns);
    for (size_t i = 0; i < numpanels; ++i) {
        const auto &p = *panels[i];
        p.Derivative_01(m_sc, tmp, tmp_n);
        Derivative_01_sc.block(0, N * i, N, N) = tmp;
        Derivative_01_sc_n.block(0, N * i, N, N) = tmp_n;
        p.Derivative_01(m_tc, tmp, tmp_n);
        Derivative_01_tc.block(0, N * i, N, N) = tmp;
        Derivative_01_tc_n.block(0, N * i, N, N) = tmp_n;
        p.Derivative_01(m_sa, tmp, tmp_n);
        Derivative_01_sa.block(0, N * i, N, N) = tmp;
        Derivative_01_sa_n.block(0, N * i, N, N) = tmp_n;
        p.Derivative_01(m_ta, tmp, tmp_n);
        Derivative_01_ta.block(0, N * i, N, N) = tmp;
        Derivative_01_ta_n.block(0, N * i, N, N) = tmp_n;
        p.Derivative_01_swapped(m_sa, tmp, tmp_n, true);
        Derivative_01_swapped_sa.block(0, N * i, N, N) = tmp;
        Derivative_01_swapped_sa_n.block(0, N * i, N, N) = tmp_n;
        p.Derivative_01_swapped(m_ta, tmp, tmp_n, true);
        Derivative_01_swapped_ta.block(0, N * i, N, N) = tmp;
        Derivative_01_swapped_ta_n.block(0, N * i, N, N) = tmp_n;
        p.Derivative_01(m_sg, tmps, tmps_n);
        Derivative_01_sg.block(0, Ns * i, Ns, Ns) = tmps;
        Derivative_01_sg_n.block(0, Ns * i, Ns, Ns) = tmps_n;
        p.Derivative_01(m_tg, tmps, tmps_n);
        Derivative_01_tg.block(0, Ns * i, Ns, Ns) = tmps;
        Derivative_01_tg_n.block(0, Ns * i, Ns, Ns) = tmps_n;
        op_sc.block(0, i * N, N, N) = p[m_sc];
        op_tc.block(0, i * N, N, N) = p[m_tc];
        op_sa.block(0, i * N, N, N) = p[m_sa];
        op_ta.block(0, i * N, N, N) = p[m_ta];
        swapped_op_sa.block(0, i * N, N, N) = p.swapped_op(m_sa);
        swapped_op_ta.block(0, i * N, N, N) = p.swapped_op(m_ta);
        op_sg.block(0, i * Ns, Ns, Ns) = p[m_sg];
        op_tg.block(0, i * Ns, Ns, Ns) = p[m_tg];
    }
}

// return true iff the panels P1 and P2 are adjacent, compute the SWAP boolean alongside
bool BuilderData::is_adjacent(size_t i, size_t j, bool &swap) const {
    const auto &p1 = *panels[i], &p2 = *panels[j];
    double t1 = (p1(1) - p2(-1)).norm() / 100., t2 = (p1(-1) - p2(1)).norm() / 100.;
    if (t1 < epsilon || t2 < epsilon) {
        swap = t1 > epsilon;
        return true;
    }
    return false;
}

GalerkinMatrixBuilder::GalerkinMatrixBuilder(const BuilderData &builder_data)
: data(builder_data)
{
    size_t N = data.getCGaussQROrder(), Ns = data.getGaussQROrder();
    size_t Qtest = data.getQtest(), Qtrial = data.getQtrial();
    size_t dim_test = data.getTestSpaceDimension(), dim_trial = data.getTrialSpaceDimension();
    m_zero.setZero(N, N);
    m_zero_s.setZero(Ns, Ns);
    for (size_t K = 0; K < 2; ++K) {
        m_h0[K].resize(N, N);
        m_h1[K].resize(N, N);
        m_v[K].resize(N, N);
        m_tangent[K].resize(N, N);
        m_tangent_p[K].resize(N, N);
        masked1[K].resize(N, N);
        masked2[K].resize(N, N);
        masked3[K].resize(N, N);
        m_v_norm[K].resize(N, N);
        m_v_norm2[K].resize(N, N);
        m_tangent_norm[K].resize(N, N);
        m_tangent_p_norm[K].resize(N, N);
    }
    m_h0_s.resize(Ns, Ns);
    m_h1_s.resize(Ns, Ns);
    m_v_s.resize(Ns, Ns);
    m_tangent_s.resize(Ns, Ns);
    m_tangent_p_s.resize(Ns, Ns);
    masked1_s.resize(Ns, Ns);
    masked2_s.resize(Ns, Ns);
    masked3_s.resize(Ns, Ns);
    m_v_norm_s.resize(Ns, Ns);
    m_v_norm2_s.resize(Ns, Ns);
    m_tangent_norm_s.resize(Ns, Ns);
    m_tangent_p_norm_s.resize(Ns, Ns);
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
    interaction_matrix_assembly_time = hankel_computation_time = 0;
}

// compute values required for the coinciding panels case
inline void GalerkinMatrixBuilder::compute_coinciding(size_t i) throw() {
    size_t I, J, N = data.getCGaussQROrder();
    m_tangent_p[0] = data.Derivative_01_tc.block(0, N * i, N, N);
    m_tangent_p_norm[0] = data.Derivative_01_tc_n.block(0, N * i, N, N);
    m_tangent[0] = data.Derivative_01_sc.block(0, N * i, N, N);
    m_tangent_norm[0] = data.Derivative_01_sc_n.block(0, N * i, N, N);
    m_v[0] = data.op_sc.block(0, N * i, N, N) - data.op_tc.block(0, N * i, N, N);
    m_v_norm2[0] = m_v[0].cwiseAbs2();
    m_v_norm[0] = m_v_norm2[0].cwiseSqrt();
    auto tic = chrono::high_resolution_clock::now();
    if (k_real_positive) {
        complex_bessel::H1_01_i(ksqrtc.real() * m_v_norm[0], m_h0[0], m_h1[0]);
    } else {
        for (J = 0; J < N; ++J) {
            for (I = 0; I < N; ++I) {
                complex_t arg = ksqrtc * m_v_norm[0](I, J);
                m_h0[0](I, J) = complex_bessel::H1(0, arg);
                m_h1[0](I, J) = complex_bessel::H1(1, arg);
            }
        }
        m_h0[0] *= 1i;
        m_h1[0] *= 1i;
    }
    auto toc = chrono::high_resolution_clock::now();
    hankel_computation_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
}

// compute values required for the adjacent panels case
inline void GalerkinMatrixBuilder::compute_adjacent(size_t i, size_t j, bool swap) throw() {
    size_t I, J, K, N = data.getCGaussQROrder();
    if (swap) {
        m_tangent_p[0] = data.Derivative_01_swapped_ta.block(0, N * j, N, N);
        m_tangent_p_norm[0] = data.Derivative_01_swapped_ta_n.block(0, N * j, N, N);
        m_tangent_p[1] = data.Derivative_01_swapped_sa.block(0, N * j, N, N);
        m_tangent_p_norm[1] = data.Derivative_01_swapped_sa_n.block(0, N * j, N, N);
        m_tangent[0] = data.Derivative_01_sa.block(0, N * i, N, N);
        m_tangent_norm[0] = data.Derivative_01_sa_n.block(0, N * i, N, N);
        m_tangent[1] = data.Derivative_01_ta.block(0, N * i, N, N);
        m_tangent_norm[1] = data.Derivative_01_ta_n.block(0, N * i, N, N);
    } else {
        m_tangent_p[0] = data.Derivative_01_ta.block(0, N * j, N, N);
        m_tangent_p_norm[0] = data.Derivative_01_ta_n.block(0, N * j, N, N);
        m_tangent_p[1] = data.Derivative_01_sa.block(0, N * j, N, N);
        m_tangent_p_norm[1] = data.Derivative_01_sa_n.block(0, N * j, N, N);
        m_tangent[0] = data.Derivative_01_swapped_sa.block(0, N * i, N, N);
        m_tangent_norm[0] = data.Derivative_01_swapped_sa_n.block(0, N * i, N, N);
        m_tangent[1] = data.Derivative_01_swapped_ta.block(0, N * i, N, N);
        m_tangent_norm[1] = data.Derivative_01_swapped_ta_n.block(0, N * i, N, N);
    }
    m_v[0] = swap ? data.op_sa.block(0, N * i, N, N) - data.swapped_op_ta.block(0, N * j, N, N)
                  : data.swapped_op_sa.block(0, N * i, N, N) - data.op_ta.block(0, N * j, N, N);
    m_v[1] = swap ? data.op_ta.block(0, N * i, N, N) - data.swapped_op_sa.block(0, N * j, N, N)
                  : data.swapped_op_ta.block(0, N * i, N, N) - data.op_sa.block(0, N * j, N, N);
    m_v_norm2[0] = m_v[0].cwiseAbs2();
    m_v_norm2[1] = m_v[1].cwiseAbs2();
    m_v_norm[0] = m_v_norm2[0].cwiseSqrt();
    m_v_norm[1] = m_v_norm2[1].cwiseSqrt();
    auto tic = chrono::high_resolution_clock::now();
    if (k_real_positive) {
        complex_bessel::H1_01_i(ksqrtc.real() * m_v_norm[0], m_h0[0], m_h1[0]);
        complex_bessel::H1_01_i(ksqrtc.real() * m_v_norm[1], m_h0[1], m_h1[1]);
    } else for (K = 0; K < 2; ++K) {
        for (J = 0; J < N; ++J) {
            for (I = 0; I < N; ++I) {
                complex_t arg = ksqrtc * m_v_norm[K](I, J);
                m_h0[K](I, J) = complex_bessel::H1(0, arg);
                m_h1[K](I, J) = complex_bessel::H1(1, arg);
            }
        }
        m_h0[K] *= 1i;
        m_h1[K] *= 1i;
    }
    auto toc = chrono::high_resolution_clock::now();
    hankel_computation_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
}

// compute values required for the disjoint panels case
inline void GalerkinMatrixBuilder::compute_general(size_t i, size_t j) throw() {
    size_t I, J, N = data.getGaussQROrder();
    m_tangent_p_s = data.Derivative_01_tg.block(0, N * j, N, N);
    m_tangent_p_norm_s = data.Derivative_01_tg_n.block(0, N * j, N, N);
    m_tangent_s = data.Derivative_01_sg.block(0, N * i, N, N);
    m_tangent_norm_s = data.Derivative_01_sg_n.block(0, N * i, N, N);
    m_v_s = data.op_sg.block(0, N * i, N, N) - data.op_tg.block(0, N * j, N, N);
    m_v_norm2_s = m_v_s.cwiseAbs2();
    m_v_norm_s = m_v_norm2_s.cwiseSqrt();
    auto tic = chrono::high_resolution_clock::now();
    if (k_real_positive) {
        complex_bessel::H1_01_i(ksqrtc.real() * m_v_norm_s, m_h0_s, m_h1_s);
    } else {
        for (J = 0; J < N; ++J) {
            for (I = 0; I < N; ++I) {
                complex_t arg = ksqrtc * m_v_norm_s(I, J);
                m_h0_s(I, J) = complex_bessel::H1(0, arg);
                m_h1_s(I, J) = complex_bessel::H1(1, arg);
            }
        }
        m_h0_s *= 1i;
        m_h1_s *= 1i;
    }
    auto toc = chrono::high_resolution_clock::now();
    hankel_computation_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
}

// compute interaction matrices for the K submatrix, coinciding panels case
void GalerkinMatrixBuilder::double_layer_coinciding(int der) throw() {
    size_t i, j, K, N = data.getCGaussQROrder(), Qtest = data.getQtest(), Qtrial = data.getQtrial();
    const auto &v = m_v[0];
    const auto &v_norm = m_v_norm[0], &v_norm2 = m_v_norm2[0];
    const auto &h1 = m_h1[0], &h0 = m_h0[0];
    auto mask = (ksqrtca * v_norm) > epsilon;
    Eigen::ArrayXXd vdotn, cf;
    double_layer_interaction_matrix.setZero();
    masked1[0] = mask.select(ksqrtc * h1 / v_norm, (v_norm > epsilon).select(M_2_PI / v_norm2, m_zero));
    if (der > 0) {
        double_layer_der_interaction_matrix.setZero();
        masked2[0] = mask.select(h0, m_zero);
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix.setZero();
        masked3[0] = masked2[0] - mask.select(h1 * ksqrtc * v_norm, m_zero);
    }
    for (K = 0; K < 2; ++K) {
        const auto &tangent = (K == 0 ? m_tangent_p : m_tangent)[0];
        vdotn = (K * 2 - 1.) * data.m_wc * (K == 0 ? m_tangent_norm : m_tangent_p_norm)[0] *
            (v.imag() * tangent.real() - v.real() * tangent.imag());
        for (j = 0; j < Qtrial; ++j) {
            for (i = 0; i < Qtest; ++i) {
                cf = vdotn * (K == 0 ? data.m_double_layer_coinciding_fg : data.m_double_layer_coinciding_fg_t).block(i * N, j * N, N, N);
                double_layer_interaction_matrix(i, j) += (cf * masked1[0]).sum();
                if (der > 0)
                    double_layer_der_interaction_matrix(i, j) += (cf * masked2[0]).sum();
                if (der > 1)
                    double_layer_der2_interaction_matrix(i, j) += (cf * masked3[0]).sum();
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
    size_t i, j, K, N = data.getCGaussQROrder(), Qtest = data.getQtest(), Qtrial = data.getQtrial();
    double_layer_interaction_matrix.setZero();
    if (der > 0)
        double_layer_der_interaction_matrix.setZero();
    if (der > 1)
        double_layer_der2_interaction_matrix.setZero();
    Eigen::ArrayXXd vdotn, cf;
    for (K = 0; K < 2; ++K) {
        const auto &tangent = (transp ? m_tangent : m_tangent_p)[K];
        const auto &v = m_v[K];
        const auto &v_norm = m_v_norm[K], &v_norm2 = m_v_norm2[K];
        const auto &h1 = m_h1[K], &h0 = m_h0[K];
        auto mask = (ksqrtca * v_norm) > epsilon;
        if (!transp) {
            if (k_real_positive)
                masked1[K] = mask.select((ksqrtc.real() / v_norm) * h1, (v_norm > epsilon).select(M_2_PI / v_norm2, m_zero));
            else masked1[K] = mask.select(ksqrtc * h1 / v_norm, (v_norm > epsilon).select(M_2_PI / v_norm2, m_zero));
            if (der > 0)
                masked2[K] = mask.select(h0, m_zero);
            if (der > 1) {
                if (k_real_positive)
                    masked3[K] = masked2[K] - mask.select((ksqrtc.real() * v_norm) * h1, m_zero);
                else masked3[K] = masked2[K] - mask.select(ksqrtc * v_norm * h1, m_zero);
            }
        }
        vdotn = (transp? 1.0 : -1.0) * data.m_wa * (v.imag() * tangent.real() - v.real() * tangent.imag()) *
            (transp ? m_tangent_p_norm : m_tangent_norm)[K];
        const Eigen::ArrayXXd &fg = swap ? ((transp ? K == 0 : K == 1) ? data.m_double_layer_adjacent_fg_swap_t
                                                                       : data.m_double_layer_adjacent_fg_swap)
                                         : ((transp ? K == 0 : K == 1) ? data.m_double_layer_adjacent_fg_t
                                                                       : data.m_double_layer_adjacent_fg);
        for (j = 0; j < Qtrial; ++j) {
            for (i = 0; i < Qtest; ++i) {
                cf = vdotn * fg.block(i * N, j * N, N, N);
                double_layer_interaction_matrix(i, j) += (cf * masked1[K]).sum();
                if (der > 0)
                    double_layer_der_interaction_matrix(i, j) += (cf * masked2[K]).sum();
                if (der > 1)
                    double_layer_der2_interaction_matrix(i, j) += (cf * masked3[K]).sum();
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
    size_t i, j, N = data.getGaussQROrder(), Qtest = data.getQtest(), Qtrial = data.getQtrial();
    Eigen::ArrayXXd cf, vdotn;
    const Eigen::ArrayXXd &fg = transp ? data.m_double_layer_general_fg_t : data.m_double_layer_general_fg;
    if (!transp) {
        vdotn = data.m_wg * (m_v_s.real() * m_tangent_p_s.imag() - m_v_s.imag() * m_tangent_p_s.real()) * m_tangent_norm_s;
        auto mask = ksqrtca * m_v_norm_s > epsilon;
        if (k_real_positive)
            masked1_s = mask.select((ksqrtc / m_v_norm_s) * m_h1_s, (m_v_norm_s > epsilon).select(M_2_PI / m_v_norm2_s, m_zero_s));
        else masked1_s = mask.select(ksqrtc * m_h1_s / m_v_norm_s, (m_v_norm_s > epsilon).select(M_2_PI / m_v_norm2_s, m_zero_s));
        if (der > 0)
            masked2_s = mask.select(m_h0_s, m_zero_s);
        if (der > 1) {
            if (k_real_positive)
                masked3_s = masked2_s - mask.select((ksqrtc.real() * m_v_norm_s) * m_h1_s, m_zero_s);
            else masked3_s = masked2_s - mask.select(ksqrtc * m_h1_s * m_v_norm_s, m_zero_s);
        }
    } else vdotn = data.m_wg * (m_v_s.imag() * m_tangent_s.real() - m_v_s.real() * m_tangent_s.imag()) * m_tangent_p_norm_s;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            cf = vdotn * fg.block(i * N, j * N, N, N);
            double_layer_interaction_matrix(i, j) = (cf * masked1_s).sum();
            if (der > 0)
                double_layer_der_interaction_matrix(i, j) = (cf * masked2_s).sum();
            if (der > 1)
                double_layer_der2_interaction_matrix(i, j) = (cf * masked3_s).sum();
        }
    }
    if (der > 0)
        double_layer_der_interaction_matrix *= k * c;
    if (der > 1)
        double_layer_der2_interaction_matrix *= c;
}

// compute interaction matrices for the W submatrix, coinciding panels case
void GalerkinMatrixBuilder::hypersingular_coinciding(int der) throw() {
    size_t i, j, K, N = data.getCGaussQROrder(), Qtrial = data.getQtrial();
    const auto &tangent_p = m_tangent_p[0], &tangent = m_tangent[0];
    const auto &h0 = m_h0[0];
    const auto &v_norm = m_v_norm[0], &v_norm2 = m_v_norm2[0];
    Eigen::ArrayXXd tdottp = 2.0 * (tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real()) * data.m_wc;
    auto mask = ksqrtca * v_norm > epsilon;
    Eigen::ArrayXXcd h1_vnorm, cf;
    Eigen::ArrayXXd temp;
    hypersingular_interaction_matrix.setZero();
    masked1[0] = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
    if (der > 0) {
        hypersingular_der_interaction_matrix.setZero();
        h1_vnorm = m_h1[0] * v_norm;
    }
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
                cf = data.m_wc * data.m_hypersingular_coinciding_fg_arc.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) -
                    kkc * (temp = data.m_hypersingular_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) * tdottp);
                hypersingular_interaction_matrix(i, j) += (cf * masked1[0]).sum();
                if (der > 0)
                    hypersingular_der_interaction_matrix(i, j) += mask.select(h1_vnorm * cf + ksqrtc * h0 * temp, m_zero).sum();
                if (der > 1)
                    hypersingular_der2_interaction_matrix(i, j) += mask.select(
                        (h0 * v_norm2 - h1_vnorm / ksqrtc) * cf - temp * (2.0 * ksqrtc * h1_vnorm - h0), m_zero).sum();
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
    size_t i, j, K, N = data.getCGaussQROrder(), Qtrial = data.getQtrial();
    hypersingular_interaction_matrix.setZero();
    if (der > 0)
        hypersingular_der_interaction_matrix.setZero();
    if (der > 1)
        hypersingular_der2_interaction_matrix.setZero();
    Eigen::ArrayXXcd h1_vnorm, cf;
    Eigen::ArrayXXd temp, tdottp;
    for (K = 0; K < 2; ++K) {
        const auto &tangent_p = m_tangent_p[K], &tangent = m_tangent[K];
        const auto &h0 = m_h0[K];
        const auto &v_norm = m_v_norm[K], &v_norm2 = m_v_norm2[K];
        tdottp = 2.0 * (tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real()) * data.m_wa;
        auto mask = ksqrtca * v_norm > epsilon;
        masked1[0] = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
        if (der > 0)
            h1_vnorm = m_h1[K] * v_norm;
        for (j = 0; j < Qtrial; ++j) {
            for (i = 0; i < Qtrial; ++i) {
                cf = data.m_wa * (swap ? (K > 0 ? data.m_hypersingular_adjacent_fg_arc.block(j * N, i * N, N, N)
                                                : data.m_hypersingular_adjacent_fg_arc_swap.block(i * N, j * N, N, N))
                                       : (K > 0 ? data.m_hypersingular_adjacent_fg_arc_swap.block(j * N, i * N, N, N)
                                                : data.m_hypersingular_adjacent_fg_arc.block(i * N, j * N, N, N))) -
                    kkc * (temp = (swap ? (K > 0 ? data.m_hypersingular_adjacent_fg.block(j * N, i * N, N, N)
                                                 : data.m_hypersingular_adjacent_fg_swap.block(i * N, j * N, N, N))
                                        : (K > 0 ? data.m_hypersingular_adjacent_fg_swap.block(j * N, i * N, N, N)
                                                 : data.m_hypersingular_adjacent_fg.block(i * N, j * N, N, N))) * tdottp);
                hypersingular_interaction_matrix(i, j) += (cf * masked1[0]).sum();
                if (der > 0)
                    hypersingular_der_interaction_matrix(i, j) += mask.select(h1_vnorm * cf + ksqrtc * h0 * temp, m_zero).sum();
                if (der > 1)
                    hypersingular_der2_interaction_matrix(i, j) += mask.select(
                        (h0 * v_norm2 - h1_vnorm / ksqrtc) * cf - temp * (2.0 * ksqrtc * h1_vnorm - h0), m_zero).sum();
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
    size_t i, j, N = data.getGaussQROrder(), Qtrial = data.getQtrial();
    const auto &tangent_p = m_tangent_p_s, &tangent = m_tangent_s;
    const auto &h0 = m_h0_s;
    const auto &v_norm = m_v_norm_s, &v_norm2 = m_v_norm2_s;
    Eigen::ArrayXXd tdottp = 2.0 * (tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real()) * data.m_wg;
    auto mask = ksqrtca * v_norm > epsilon;
    Eigen::ArrayXXcd h1_vnorm, cf;
    Eigen::ArrayXXd temp;
    masked1_s = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero_s));
    if (der > 0)
        h1_vnorm = m_h1_s * v_norm;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            cf = data.m_wg * data.m_hypersingular_general_fg_arc.block(i * N, j * N, N, N) -
                kkc * (temp = data.m_hypersingular_general_fg.block(i * N, j * N, N, N) * tdottp);
            hypersingular_interaction_matrix(i, j) = (cf * masked1_s).sum();
            if (der > 0)
                hypersingular_der_interaction_matrix(i, j) =
                    mask.select(h1_vnorm * cf + ksqrtc * h0 * temp, m_zero_s).sum();
            if (der > 1)
                hypersingular_der2_interaction_matrix(i, j) = mask.select(
                    (h0 * v_norm2 - h1_vnorm / ksqrtc) * cf - temp * (2.0 * ksqrtc * h1_vnorm - h0), m_zero_s).sum();
        }
    }
    if (der > 0)
        hypersingular_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        hypersingular_der2_interaction_matrix *= -c;
}

// compute interaction matrices for the V submatrix, coinciding panels case
void GalerkinMatrixBuilder::single_layer_coinciding(int der) throw() {
    size_t i, j, K, N = data.getCGaussQROrder(), Qtest = data.getQtest();
    const auto &v_norm = m_v_norm[0];
    const auto &h0 = m_h0[0], &h1 = m_h1[0];
    Eigen::ArrayXXd cf, ttp_norm = m_tangent_norm[0] * m_tangent_p_norm[0] * data.m_wc;
    auto mask = ksqrtca * v_norm > epsilon;
    single_layer_interaction_matrix.setZero();
    masked1[0] = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
    if (der > 0) {
        single_layer_der_interaction_matrix.setZero();
        masked2[0] = mask.select(h1 * v_norm, m_zero);
    }
    if (der > 1) {
        single_layer_der2_interaction_matrix.setZero();
        masked3[0] = mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm, m_zero);
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
                cf = ttp_norm * data.m_single_layer_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N);
                single_layer_interaction_matrix(i, j) += (cf * masked1[0]).sum();
                if (der > 0)
                    single_layer_der_interaction_matrix(i, j) += (cf * masked2[0]).sum();
                if (der > 1)
                    single_layer_der2_interaction_matrix(i, j) += (cf * masked3[0]).sum();
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
    size_t i, j, K, N = data.getCGaussQROrder(), Qtest = data.getQtest();
    single_layer_interaction_matrix.setZero();
    if (der > 0)
        single_layer_der_interaction_matrix.setZero();
    if (der > 1)
        single_layer_der2_interaction_matrix.setZero();
    Eigen::ArrayXXd ttp_norm, cf;
    for (K = 0; K < 2; ++K) {
        const auto &v_norm = m_v_norm[K];
        const auto &h0 = m_h0[K], &h1 = m_h1[K];
        auto mask = ksqrtca * v_norm > epsilon;
        masked1[K] = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero));
        if (der > 0)
            masked2[K] = mask.select(h1 * v_norm, m_zero);
        if (der > 1)
            masked3[K] = mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm, m_zero);
        ttp_norm = m_tangent_norm[K] * m_tangent_p_norm[K] * data.m_wa;
        for (j = 0; j < Qtest; ++j) {
            for (i = 0; i < Qtest; ++i) {
                cf = ttp_norm * (swap ? (K > 0 ? data.m_single_layer_adjacent_fg.block(j * N, i * N, N, N)
                                               : data.m_single_layer_adjacent_fg_swap.block(i * N, j * N, N, N))
                                      : (K > 0 ? data.m_single_layer_adjacent_fg_swap.block(j * N, i * N, N, N)
                                               : data.m_single_layer_adjacent_fg.block(i * N, j * N, N, N)));
                single_layer_interaction_matrix(i, j) += (cf * masked1[K]).sum();
                if (der > 0)
                    single_layer_der_interaction_matrix(i, j) += (cf * masked2[K]).sum();
                if (der > 1)
                    single_layer_der2_interaction_matrix(i, j) += (cf * masked3[K]).sum();
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
    size_t i, j, N = data.getGaussQROrder(), Qtest = data.getQtest();
    const auto &v_norm = m_v_norm_s;
    const auto &h0 = m_h0_s, &h1 = m_h1_s;
    Eigen::ArrayXXd cf, ttp_norm = m_tangent_norm_s * m_tangent_p_norm_s * data.m_wg;
    auto mask = ksqrtca * v_norm > epsilon;
    masked1_s = mask.select(h0, (v_norm > epsilon).select(-M_2_PI * v_norm.log(), m_zero_s));
    if (der > 0)
        masked2_s = mask.select(h1 * v_norm, m_zero_s);
    if (der > 1)
        masked3_s = mask.select((h1 / ksqrtc - h0 * v_norm) * v_norm, m_zero_s);
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            cf = ttp_norm * data.m_single_layer_general_fg.block(i * N, j * N, N, N);
            single_layer_interaction_matrix(i, j) = (cf * masked1_s).sum();
            if (der > 0)
                single_layer_der_interaction_matrix(i, j) = (cf * masked2_s).sum();
            if (der > 1)
                single_layer_der2_interaction_matrix(i, j) = (cf * masked3_s).sum();
        }
    }
    if (der > 0)
        single_layer_der_interaction_matrix *= -sqrtc;
    if (der > 1)
        single_layer_der2_interaction_matrix *= c;
}

void GalerkinMatrixBuilder::all_coinciding(int der) throw() {
    size_t i, j, K, N = data.getCGaussQROrder(), Q = data.getQtest();
    const auto &v = m_v[0];
    const auto &v_norm = m_v_norm[0], &v_norm2 = m_v_norm2[0];
    const auto &h1 = m_h1[0], &h0 = m_h0[0];
    const auto &tangent_p = m_tangent_p[0], &tangent = m_tangent[0];
    Eigen::ArrayXXd tdottp = 2.0 * (tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real()) * data.m_wc;
    Eigen::ArrayXXd ttp_norm = m_tangent_norm[0] * m_tangent_p_norm[0] * data.m_wc;
    auto mask = (ksqrtca * v_norm) > epsilon;
    auto mask2 = v_norm > epsilon;
    masked1[0] = mask.select(ksqrtc * h1 / v_norm, mask2.select(M_2_PI / v_norm2, m_zero));
    Eigen::ArrayXXcd masked1_hg = mask.select(h0, mask2.select(-M_2_PI * v_norm.log(), m_zero));
    Eigen::ArrayXXcd masked2_g, masked3_g, cf;
    Eigen::ArrayXXd temp, vdotn, cfr1, cfr2;
    double_layer_interaction_matrix.setZero();
    hypersingular_interaction_matrix.setZero();
    single_layer_interaction_matrix.setZero();
    if (der > 0) {
        double_layer_der_interaction_matrix.setZero();
        hypersingular_der_interaction_matrix.setZero();
        single_layer_der_interaction_matrix.setZero();
        masked2[0] = mask.select(h0, m_zero);
        masked2_g = mask.select(h1 * v_norm, m_zero);
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix.setZero();
        hypersingular_der2_interaction_matrix.setZero();
        single_layer_der2_interaction_matrix.setZero();
        masked3[0] = masked2[0] - ksqrtc * masked2_g;
        masked3_g = masked2_g / ksqrtc - masked2[0] * v_norm2;
    }
    for (K = 0; K < 2; ++K) {
        vdotn = (K * 2 - 1.) * data.m_wc * (K == 0 ? m_tangent_norm : m_tangent_p_norm)[0] *
            (v.imag() * (K == 0 ? tangent_p : tangent).real() - v.real() * (K == 0 ? tangent_p : tangent).imag());
        const Eigen::ArrayXXd &fg = K == 0 ? data.m_double_layer_coinciding_fg : data.m_double_layer_coinciding_fg_t;
        for (j = 0; j < Q; ++j) {
            for (i = 0; i < Q; ++i) {
                cfr1 = vdotn * fg.block(i * N, j * N, N, N);
                double_layer_interaction_matrix(i, j) += (cfr1 * masked1[0]).sum();
                if (der > 0)
                    double_layer_der_interaction_matrix(i, j) += (cfr1 * masked2[0]).sum();
                if (der > 1)
                    double_layer_der2_interaction_matrix(i, j) += (cfr1 * masked3[0]).sum();
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
                    temp = data.m_hypersingular_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) * tdottp;
                    if (k_real_positive)
                        cfr1 = data.m_wc * data.m_hypersingular_coinciding_fg_arc.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) - kkc.real() * temp;
                    else cf = data.m_wc * data.m_hypersingular_coinciding_fg_arc.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N) - kkc * temp;
                    cfr2 = ttp_norm * data.m_single_layer_coinciding_fg.block((K > 0 ? j : i) * N, (K > 0 ? i : j) * N, N, N);
                    single_layer_interaction_matrix(i, j) += (cfr2 * masked1_hg).sum();
                    if (k_real_positive)
                        hypersingular_interaction_matrix(i, j) += (cfr1 * masked1_hg).sum();
                    else hypersingular_interaction_matrix(i, j) += (cf * masked1_hg).sum();
                    if (der > 0) {
                        single_layer_der_interaction_matrix(i, j) += (cfr2 * masked2_g).sum();
                        if (k_real_positive)
                            hypersingular_der_interaction_matrix(i, j) += (masked2_g * cfr1 + (ksqrtc.real() * temp) * masked2[0]).sum();
                        else hypersingular_der_interaction_matrix(i, j) += (masked2_g * cf + ksqrtc * masked2[0] * temp).sum();
                    }
                    if (der > 1) {
                        single_layer_der2_interaction_matrix(i, j) += (cfr2 * masked3_g).sum();
                        if (k_real_positive)
                            hypersingular_der2_interaction_matrix(i, j) -= (cfr1 * masked3_g + temp * (masked2[0] - 2.0 * masked3[0])).sum();
                        else hypersingular_der2_interaction_matrix(i, j) -= (masked3_g * cf + temp * (masked2[0] - 2.0 * masked3[0])).sum();
                    }
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

void GalerkinMatrixBuilder::all_adjacent(bool swap, int der) throw() {
    size_t i, j, K, N = data.getCGaussQROrder(), Q = data.getQtest();
    double_layer_interaction_matrix.setZero();
    hypersingular_interaction_matrix.setZero();
    single_layer_interaction_matrix.setZero();
    if (der > 0) {
        double_layer_der_interaction_matrix.setZero();
        hypersingular_der_interaction_matrix.setZero();
        single_layer_der_interaction_matrix.setZero();
    }
    if (der > 1) {
        double_layer_der2_interaction_matrix.setZero();
        hypersingular_der2_interaction_matrix.setZero();
        single_layer_der2_interaction_matrix.setZero();
    }
    Eigen::ArrayXXcd masked2_g, masked3_g, cf;
    Eigen::ArrayXXd temp, cfr1, cfr2, cfr3;
    for (K = 0; K < 2; ++K) {
        const auto &tangent_p = m_tangent_p[K], &tangent = m_tangent[K];
        const auto &v = m_v[K];
        const auto &v_norm = m_v_norm[K], &v_norm2 = m_v_norm2[K];
        const auto &h1 = m_h1[K], &h0 = m_h0[K];
        Eigen::ArrayXXd tdottp = 2.0 * (tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real()) * data.m_wa;
        Eigen::ArrayXXd vdotn = data.m_wa * (v.real() * tangent_p.imag() - v.imag() * tangent_p.real()) * m_tangent_norm[K];
        Eigen::ArrayXXd ttp_norm = m_tangent_norm[K] * m_tangent_p_norm[K] * data.m_wa;
        auto mask = (ksqrtca * v_norm) > epsilon;
        auto mask2 = v_norm > epsilon;
        Eigen::ArrayXXcd masked1_hg = mask.select(h0, mask2.select(-M_2_PI * v_norm.log(), m_zero));
        masked1[K] = mask.select(ksqrtc * h1 / v_norm, mask2.select(M_2_PI / v_norm2, m_zero));
        if (der > 0) {
            masked2[K] = mask.select(h0, m_zero);
            masked2_g = mask.select(h1 * v_norm, m_zero);
        }
        if (der > 1) {
            masked3[K] = masked2[K] - ksqrtc * masked2_g;
            masked3_g = masked2_g / ksqrtc - masked2[K] * v_norm2;
        }
        const Eigen::ArrayXXd &fg = swap ? (K == 1 ? data.m_double_layer_adjacent_fg_swap_t : data.m_double_layer_adjacent_fg_swap)
                                         : (K == 1 ? data.m_double_layer_adjacent_fg_t : data.m_double_layer_adjacent_fg);
        for (j = 0; j < Q; ++j) {
            for (i = 0; i < Q; ++i) {
                cfr1 = vdotn * fg.block(i * N, j * N, N, N);
                cfr2 = ttp_norm * (swap ? (K > 0 ? data.m_single_layer_adjacent_fg.block(j * N, i * N, N, N)
                                                 : data.m_single_layer_adjacent_fg_swap.block(i * N, j * N, N, N))
                                        : (K > 0 ? data.m_single_layer_adjacent_fg_swap.block(j * N, i * N, N, N)
                                                 : data.m_single_layer_adjacent_fg.block(i * N, j * N, N, N)));
                temp = (swap ? (K > 0 ? data.m_hypersingular_adjacent_fg.block(j * N, i * N, N, N)
                                      : data.m_hypersingular_adjacent_fg_swap.block(i * N, j * N, N, N))
                             : (K > 0 ? data.m_hypersingular_adjacent_fg_swap.block(j * N, i * N, N, N)
                                      : data.m_hypersingular_adjacent_fg.block(i * N, j * N, N, N))) * tdottp;
                if (k_real_positive)
                    cfr3 = data.m_wa * (swap ? (K > 0 ? data.m_hypersingular_adjacent_fg_arc.block(j * N, i * N, N, N)
                                                      : data.m_hypersingular_adjacent_fg_arc_swap.block(i * N, j * N, N, N))
                                             : (K > 0 ? data.m_hypersingular_adjacent_fg_arc_swap.block(j * N, i * N, N, N)
                                                      : data.m_hypersingular_adjacent_fg_arc.block(i * N, j * N, N, N))) - kkc.real() * temp;
                else cf = data.m_wa * (swap ? (K > 0 ? data.m_hypersingular_adjacent_fg_arc.block(j * N, i * N, N, N)
                                                     : data.m_hypersingular_adjacent_fg_arc_swap.block(i * N, j * N, N, N))
                                            : (K > 0 ? data.m_hypersingular_adjacent_fg_arc_swap.block(j * N, i * N, N, N)
                                                     : data.m_hypersingular_adjacent_fg_arc.block(i * N, j * N, N, N))) - kkc * temp;
                double_layer_interaction_matrix(i, j) += (cfr1 * masked1[K]).sum();
                single_layer_interaction_matrix(i, j) += (cfr2 * masked1_hg).sum();
                if (k_real_positive)
                    hypersingular_interaction_matrix(i, j) += (cfr3 * masked1_hg).sum();
                else hypersingular_interaction_matrix(i, j) += (cf * masked1_hg).sum();
                if (der > 0) {
                    double_layer_der_interaction_matrix(i, j) += (cfr1 * masked2[K]).sum();
                    single_layer_der_interaction_matrix(i, j) += (cfr2 * masked2_g).sum();
                    if (k_real_positive)
                        hypersingular_der_interaction_matrix(i, j) += (masked2_g * cfr3 + (ksqrtc.real() * temp) * masked2[K]).sum();
                    else hypersingular_der_interaction_matrix(i, j) += (masked2_g * cf + ksqrtc * masked2[K] * temp).sum();
                }
                if (der > 1) {
                    double_layer_der2_interaction_matrix(i, j) += (cfr1 * masked3[K]).sum();
                    single_layer_der2_interaction_matrix(i, j) += (cfr2 * masked3_g).sum();
                    if (k_real_positive)
                        hypersingular_der2_interaction_matrix(i, j) -= (masked3_g * cfr3 + temp * (masked2[K] - 2.0 * masked3[K])).sum();
                    else hypersingular_der2_interaction_matrix(i, j) -= (masked3_g * cf + temp * (masked2[K] - 2.0 * masked3[K])).sum();
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
    size_t i, j, N = data.getGaussQROrder(), Q = data.getQtest();
    Eigen::ArrayXXd ttp_norm = m_tangent_norm_s * m_tangent_p_norm_s * data.m_wg;
    Eigen::ArrayXXd tdottp = 2.0 * (m_tangent_s.imag() * m_tangent_p_s.imag() + m_tangent_s.real() * m_tangent_p_s.real()) * data.m_wg;
    auto mask = ksqrtca * m_v_norm_s > epsilon;
    auto mask2 = m_v_norm_s > epsilon;
    masked1_s = mask.select(ksqrtc * m_h1_s / m_v_norm_s, mask2.select(M_2_PI / m_v_norm2_s, m_zero_s));
    Eigen::ArrayXXcd masked1_hg = mask.select(m_h0_s, mask2.select(-M_2_PI * m_v_norm_s.log(), m_zero_s));
    Eigen::ArrayXXd vdotn = data.m_wg * (m_v_s.real() * m_tangent_p_s.imag() - m_v_s.imag() * m_tangent_p_s.real()) * m_tangent_norm_s;
    Eigen::ArrayXXcd masked2_g, masked3_g, cf;
    Eigen::ArrayXXd temp, cfr1, cfr2, cfr3;
    if (der > 0) {
        masked2_s = mask.select(m_h0_s, m_zero_s);
        masked2_g = mask.select(m_h1_s * m_v_norm_s, m_zero_s);
    }
    if (der > 1) {
        masked3_s = masked2_s - ksqrtc * masked2_g;
        masked3_g = masked2_g / ksqrtc - masked2_s * m_v_norm2_s;
    }
    for (j = 0; j < Q; ++j) {
        for (i = 0; i < Q; ++i) {
            cfr1 = vdotn * data.m_double_layer_general_fg.block(i * N, j * N, N, N);
            cfr2 = ttp_norm * data.m_single_layer_general_fg.block(i * N, j * N, N, N);
            temp = data.m_hypersingular_general_fg.block(i * N, j * N, N, N) * tdottp;
            if (k_real_positive)
                cfr3 = data.m_wg * data.m_hypersingular_general_fg_arc.block(i * N, j * N, N, N) - kkc.real() * temp;
            else cf = data.m_wg * data.m_hypersingular_general_fg_arc.block(i * N, j * N, N, N) - kkc * temp;
            double_layer_interaction_matrix(i, j) = (cfr1 * masked1_s).sum();
            single_layer_interaction_matrix(i, j) = (cfr2 * masked1_hg).sum();
            if (k_real_positive)
                hypersingular_interaction_matrix(i, j) = (cfr3 * masked1_hg).sum();
            else hypersingular_interaction_matrix(i, j) = (cf * masked1_hg).sum();
            if (der > 0) {
                double_layer_der_interaction_matrix(i, j) = (cfr1 * masked2_s).sum();
                single_layer_der_interaction_matrix(i, j) = (cfr2 * masked2_g).sum();
                if (k_real_positive)
                    hypersingular_der_interaction_matrix(i, j) = (cfr3 * masked2_g + (ksqrtc.real() * temp) * masked2_s).sum();
                else hypersingular_der_interaction_matrix(i, j) = (masked2_g * cf + ksqrtc * masked2_s * temp).sum();
            }
            if (der > 1) {
                double_layer_der2_interaction_matrix(i, j) = (cfr1 * masked3_s).sum();
                single_layer_der2_interaction_matrix(i, j) = (cfr2 * masked3_g).sum();
                if (k_real_positive)
                    hypersingular_der2_interaction_matrix(i, j) = -(cfr3 * masked3_g + temp * (masked2_s - 2.0 * masked3_s)).sum();
                else hypersingular_der2_interaction_matrix(i, j) = -(masked3_g * cf + temp * (masked2_s - 2.0 * masked3_s)).sum();
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

// initialize wavenumber and refraction index with related values
void GalerkinMatrixBuilder::initialize_parameters(const std::complex<double>& k_in, double c_in) {
    if (c_in < 1.)
        throw std::runtime_error("Refraction index must not be smaller than 1!");
    k = k_in;
    c = c_in;
    sqrtc = sqrt(c);
    ksqrtc = k * sqrtc;
    ksqrtca = std::abs(ksqrtc);
    kkc = ksqrtc * ksqrtc * .5;
    k_real_positive = k.imag() == 0. && k.real() > 0;
}

void GalerkinMatrixBuilder::assembleDoubleLayer(const std::complex<double>& k_in, double c_in, int der) {
    initialize_parameters(k_in, c_in);
    double_layer_matrix.setZero();
    if (der > 0)
        double_layer_der_matrix.setZero();
    if (der > 1)
        double_layer_der2_matrix.setZero();
    size_t i, j, I, J, II, JJ;
    size_t dim_test = data.getTestSpaceDimension(), dim_trial = data.getTrialSpaceDimension(), Qtest = data.getQtest(), Qtrial = data.getQtrial();
    bool swap;
    // Panel oriented assembly
    for (i = 0; i < dim_test; ++i) {
        for (j = 0; j < dim_trial; ++j) {
            auto tic = chrono::high_resolution_clock::now();
            if (i == j) {
                // coinciding panels
                compute_coinciding(i);
                double_layer_coinciding(der);
            } else if (data.is_adjacent(i, j, swap)) {
                // adjacent panels
                compute_adjacent(i, j, swap);
                double_layer_adjacent(swap, der, false);
            } else {
                // disjoint panels
                compute_general(i, j);
                double_layer_general(der, false);
            }
            auto toc = chrono::high_resolution_clock::now();
            interaction_matrix_assembly_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
            // Local to global mapping of the elements in interaction matrix
            for (I = 0; I < Qtest; ++I) {
                for (J = 0; J < Qtrial; ++J) {
                    II = data.test_space.LocGlobMap(I + 1, i + 1, dim_test) - 1;
                    JJ = data.trial_space.LocGlobMap(J + 1, j + 1, dim_trial) - 1;
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
    size_t dim_trial = data.getTrialSpaceDimension(), Qtrial = data.getQtrial();
    bool swap;
    // Panel oriented assembly
    for (i = 0; i < dim_trial; ++i) {
        for (j = 0; j <= i; ++j) {
            auto tic = chrono::high_resolution_clock::now();
            if (i == j) {
                // coinciding panels
                compute_coinciding(i);
                hypersingular_coinciding(der);
            } else if (data.is_adjacent(i, j, swap)) {
                // adjacent panels
                compute_adjacent(i, j, swap);
                hypersingular_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(i, j);
                hypersingular_general(der);
            }
            auto toc = chrono::high_resolution_clock::now();
            interaction_matrix_assembly_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
            // Local to global mapping of the elements in interaction matrix
            for (I = 0; I < Qtrial; ++I) {
                for (J = 0; J < Qtrial; ++J) {
                    II = data.trial_space.LocGlobMap(I + 1, i + 1, dim_trial) - 1;
                    JJ = data.trial_space.LocGlobMap(J + 1, j + 1, dim_trial) - 1;
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
                    II = data.trial_space.LocGlobMap(I + 1, j + 1, dim_trial) - 1;
                    JJ = data.trial_space.LocGlobMap(J + 1, i + 1, dim_trial) - 1;
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
    size_t dim_test = data.getTestSpaceDimension(), Qtest = data.getQtest();
    bool swap;
    // Panel oriented assembly
    for (i = 0; i < dim_test; ++i) {
        for (j = 0; j <= i; ++j) {
            auto tic = chrono::high_resolution_clock::now();
            if (i == j) {
                // coinciding panels
                compute_coinciding(i);
                single_layer_coinciding(der);
            } else if (data.is_adjacent(i, j, swap)) {
                // adjacent panels
                compute_adjacent(i, j, swap);
                single_layer_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(i, j);
                single_layer_general(der);
            }
            auto toc = chrono::high_resolution_clock::now();
            interaction_matrix_assembly_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
            // Local to global mapping of the elements in interaction matrix
            for (I = 0; I < Qtest; ++I) {
                for (J = 0; J < Qtest; ++J) {
                    II = data.test_space.LocGlobMap(I + 1, i + 1, dim_test) - 1;
                    JJ = data.test_space.LocGlobMap(J + 1, j + 1, dim_test) - 1;
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
                    II = data.test_space.LocGlobMap(I + 1, j + 1, dim_test) - 1;
                    JJ = data.test_space.LocGlobMap(J + 1, i + 1, dim_test) - 1;
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
    if (!data.testTrialSpacesAreEqual())
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
    size_t i, j, I, J, II, JJ;
    size_t dim = data.getTestSpaceDimension(), Q = data.getQtest();
    bool swap, adj;
    // Panel oriented assembly
    for (i = 0; i < dim; ++i) {
        for (j = 0; j <= i; ++j) {
            auto tic = chrono::high_resolution_clock::now();
            if (i == j) {
                // coinciding panels
                compute_coinciding(i);
                all_coinciding(der);
            } else if ((adj = data.is_adjacent(i, j, swap))) {
                // adjacent panels
                compute_adjacent(i, j, swap);
                all_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(i, j);
                all_general(der);
            }
            auto toc = chrono::high_resolution_clock::now();
            interaction_matrix_assembly_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
            // Local to global mapping of the elements in interaction matrix
            for (I = 0; I < Q; ++I) {
                for (J = 0; J < Q; ++J) {
                    II = data.test_space.LocGlobMap(I + 1, i + 1, dim) - 1;
                    JJ = data.test_space.LocGlobMap(J + 1, j + 1, dim) - 1;
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
            tic = chrono::high_resolution_clock::now();
            if (adj)
                double_layer_adjacent(!swap, der, true);
            else
                double_layer_general(der, true);
            toc = chrono::high_resolution_clock::now();
            interaction_matrix_assembly_time += chrono::duration_cast<chrono::microseconds>(toc - tic).count();
            for (I = 0; I < Q; ++I) {
                for (J = 0; J < Q; ++J) {
                    II = data.test_space.LocGlobMap(I + 1, j + 1, dim) - 1;
                    JJ = data.test_space.LocGlobMap(J + 1, i + 1, dim) - 1;
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

unsigned int GalerkinMatrixBuilder::getInteractionMatrixAssemblyTime() {
    auto ret = interaction_matrix_assembly_time;
    interaction_matrix_assembly_time = 0;
    return ret;
}

unsigned int GalerkinMatrixBuilder::getHankelComputationTime() {
    auto ret = hankel_computation_time;
    hankel_computation_time = 0;
    return ret;
}
