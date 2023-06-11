#include "galerkin_matrix_builder.hpp"
#include "cbessel.hpp"
#include <numeric>
#include <execution>

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
    numpanels = mesh.getNumPanels();
    dim_test = test_space.getSpaceDim(numpanels);
    dim_trial = trial_space_in.getSpaceDim(numpanels);
    Qtest = test_space.getQ();
    Qtrial = trial_space.getQ();
    size_t N = CGaussQR.n;
    m_h0.resize(N, 2 * N);
    m_h1.resize(N, 2 * N);
    m_v.resize(N, 2 * N);
    m_v_norm.resize(N, 2 * N);
    m_tangent.resize(N, 2 * N);
    m_tangent_p.resize(N, 2 * N);
    m_tangent_norm.resize(N, 2 * N);
    m_tangent_p_norm.resize(N, 2 * N);
    indN2.resize(2 * N);
    std::iota(indN2.begin(), indN2.end(), 0);
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
inline void GalerkinMatrixBuilder::compute_coinciding(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p) throw() {
    size_t N = CGaussQR.n;
    std::for_each(std::execution::par_unseq, indN2.cbegin(), indN2.cend(), [&](size_t J) {
        for (size_t I = 0; I < N; ++I) {
            double t = CGaussQR.x(J % N), s = t * (1. - CGaussQR.x(I));
            if (J >= N)
                std::swap(s, t);
            Eigen::Vector2d tangent_p = pi_p.Derivative_01(t);
            Eigen::Vector2d tangent = pi.Derivative_01(s);
            Eigen::Vector2d v = pi[s] - pi_p[t];
            double v_norm = v.norm();
            complex_t arg = ksqrtc * v_norm;
            m_tangent_p(I, J) = complex_t(tangent_p(0), tangent_p(1));
            m_tangent(I, J) = complex_t(tangent(0), tangent(1));
            m_tangent_p_norm(I, J) = tangent_p.norm();
            m_tangent_norm(I, J) = tangent.norm();
            m_v(I, J) = complex_t(v(0), v(1));
            m_v_norm(I, J) = v_norm;
            m_h0(I, J) = complex_bessel::H1_0_i(arg);
            m_h1(I, J) = complex_bessel::H1_1_i(arg);
        }
    });
}

// compute values required for the adjacent panels case
inline void GalerkinMatrixBuilder::compute_adjacent(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p, bool swap) throw() {
    size_t N = CGaussQR.n;
    std::for_each(std::execution::par_unseq, indN2.cbegin(), indN2.cend(), [&](size_t J) {
        for (size_t I = 0; I < N; ++I) {
            double t = CGaussQR.x(I), s = t * CGaussQR.x(J % N);
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
            m_h0(I, J) = complex_bessel::H1_0_i(arg);
            m_h1(I, J) = complex_bessel::H1_1_i(arg);
        }
    });
}

// compute values required for the disjoint panels case
inline void GalerkinMatrixBuilder::compute_general(const AbstractParametrizedCurve &pi, const AbstractParametrizedCurve &pi_p) throw() {
    size_t I, J, N = GaussQR.n;
    auto tangent_p_block = m_tangent_p.block(0, 0, N, N);
    auto tangent_block = m_tangent.block(0, 0, N, N);
    auto tangent_p_norm_block = m_tangent_p_norm.block(0, 0, N, N);
    auto tangent_norm_block = m_tangent_norm.block(0, 0, N, N);
    auto v_block = m_v.block(0, 0, N, N);
    auto v_norm_block = m_v_norm.block(0, 0, N, N);
    auto h0_block = m_h0.block(0, 0, N, N);
    auto h1_block = m_h1.block(0, 0, N, N);
    for (J = 0; J < N; ++J) {
        for (I = 0; I < N; ++I) {
            double t = GaussQR.x(J), s = GaussQR.x(I);
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
            h0_block(I, J) = complex_bessel::H1_0_i(arg);
            h1_block(I, J) = complex_bessel::H1_1_i(arg);
        }
    }
}

// compute interaction matrices for the K submatrix, coinciding panels case
inline void GalerkinMatrixBuilder::double_layer_coinciding(int der) throw() {
    size_t i, j, I, J, N = CGaussQR.n;
    double s, t, w, fg, vdotn, v_norm, cf;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < 2 * N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = CGaussQR.x(J % N), s = t * (1. - CGaussQR.x(I)), w = t * CGaussQR.w(I) * CGaussQR.w(J % N);
                    if (J >= N)
                        std::swap(s, t);
                    fg = test_space.evaluateShapeFunction(i, s) * trial_space.evaluateShapeFunction(j, t);
                    const complex_t &tangent_p = m_tangent_p(I, J), &v = m_v(I, J), &h1 = m_h1(I, J), &h0 = m_h0(I, J);
                    vdotn = v.real() * tangent_p.imag() - v.imag() * tangent_p.real();
                    cf = w * m_tangent_norm(I, J) * vdotn * fg, v_norm = m_v_norm(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += ksqrtc * h1 * cf / v_norm;
                        if (der > 0)
                            integral_der += ksqrtc * h0 * cf;
                        if (der > 1)
                            integral_der2 += (h0 - h1 * ksqrtc * v_norm) * cf;
                    } else if (v_norm > epsilon)
                        integral += M_2_PI * cf / (v_norm * v_norm);
                }
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
inline void GalerkinMatrixBuilder::double_layer_adjacent(bool swap, int der) throw() {
    size_t i, j, I, J, N = CGaussQR.n;
    double s, t, w, fg, vdotn, v_norm, cf;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < 2 * N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = CGaussQR.x(I), s = t * CGaussQR.x(J % N), w = t * CGaussQR.w(I) * CGaussQR.w(J % N);
                    if (J >= N)
                        std::swap(s, t);
                    fg = swap ? trial_space.evaluateShapeFunction_01_swapped(j, t) * test_space.evaluateShapeFunction(i, s)
                              : trial_space.evaluateShapeFunction(j, t) * test_space.evaluateShapeFunction_01_swapped(i, s);
                    const complex_t &tangent_p = m_tangent_p(I, J), &v = m_v(I, J), &h1 = m_h1(I, J), &h0 = m_h0(I, J);
                    vdotn = v.real() * tangent_p.imag() - v.imag() * tangent_p.real();
                    cf = w * m_tangent_norm(I, J) * vdotn * fg, v_norm = m_v_norm(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += ksqrtc * h1 * cf / v_norm;
                        if (der > 0)
                            integral_der += ksqrtc * h0 * cf;
                        if (der > 1)
                            integral_der2 += (h0 - h1 * ksqrtc * v_norm) * cf;
                    } else if (v_norm > epsilon)
                        integral += M_2_PI * cf / (v_norm * v_norm);
                }
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
inline void GalerkinMatrixBuilder::double_layer_general(int der) throw() {
    size_t i, j, I, J, N = GaussQR.n;
    double s, t, w, fg, vdotn, v_norm, cf;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = GaussQR.x(J), s = GaussQR.x(I), w = GaussQR.w(I) * GaussQR.w(J);
                    fg = test_space.evaluateShapeFunction(i, s) * trial_space.evaluateShapeFunction(j, t);
                    const complex_t &tangent_p = m_tangent_p(I, J), &v = m_v(I, J), &h1 = m_h1(I, J), &h0 = m_h0(I, J);
                    vdotn = v.real() * tangent_p.imag() - v.imag() * tangent_p.real();
                    cf = w * m_tangent_norm(I, J) * vdotn * fg, v_norm = m_v_norm(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += ksqrtc * h1 * cf / v_norm;
                        if (der > 0)
                            integral_der += ksqrtc * h0 * cf;
                        if (der > 1)
                            integral_der2 += (h0 - h1 * ksqrtc * v_norm) * cf;
                    } else if (v_norm > epsilon)
                        integral += M_2_PI * cf / (v_norm * v_norm);
                }
            }
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
    size_t i, j, I, J, N = CGaussQR.n;
    double t, s, w, fg, fg_arc, v_norm;
    complex_t cf;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < 2 * N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = CGaussQR.x(J % N), s = t * (1. - CGaussQR.x(I)), w = t * CGaussQR.w(I) * CGaussQR.w(J % N);
                    if (J >= N)
                        std::swap(s, t);
                    fg = trial_space.evaluateShapeFunction(i, s) * trial_space.evaluateShapeFunction(j, t);
                    fg_arc = trial_space.evaluateShapeFunctionDot_01(i, s) * trial_space.evaluateShapeFunctionDot_01(j, t);
                    const complex_t &tangent_p = m_tangent_p(I, J), &tangent = m_tangent(I, J), &h0 = m_h0(I, J), &h1 = m_h1(I, J);
                    fg *= tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
                    cf = fg_arc - kkc * fg, v_norm = m_v_norm(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += w * h0 * cf;
                        if (der > 0)
                            integral_der += w * (h1 * v_norm * cf + 2.0 * h0 * ksqrtc * fg);
                        if (der > 1)
                            integral_der2 += w * ((h0 * v_norm - h1) * cf - 2.0 * fg * (2.0 * h1 * ksqrtc - h0)) * v_norm;
                    } else if (v_norm > epsilon)
                        integral -= w * M_2_PI * log(v_norm) * cf;
                }
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
    size_t i, j, I, J, N = CGaussQR.n;
    double t, s, w, fg, fg_arc, v_norm;
    complex_t cf;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < 2 * N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = CGaussQR.x(I), s = t * CGaussQR.x(J % N), w = t * CGaussQR.w(I) * CGaussQR.w(J % N);
                    if (J >= N)
                        std::swap(s, t);
                    fg = swap ? trial_space.evaluateShapeFunction(i, s) * trial_space.evaluateShapeFunction_01_swapped(j, t)
                              : trial_space.evaluateShapeFunction_01_swapped(i, s) * trial_space.evaluateShapeFunction(j, t);
                    fg_arc = swap ? trial_space.evaluateShapeFunctionDot_01(i, s) * trial_space.evaluateShapeFunctionDot_01_swapped(j, t)
                                  : trial_space.evaluateShapeFunctionDot_01_swapped(i, s) * trial_space.evaluateShapeFunctionDot_01(j, t);
                    const complex_t &tangent_p = m_tangent_p(I, J), &tangent = m_tangent(I, J), &h0 = m_h0(I, J), &h1 = m_h1(I, J);
                    fg *= tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
                    cf = fg_arc - kkc * fg, v_norm = m_v_norm(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += w * h0 * cf;
                        if (der > 0)
                            integral_der += w * (h1 * v_norm * cf + 2.0 * h0 * ksqrtc * fg);
                        if (der > 1)
                            integral_der2 += w * ((h0 * v_norm - h1) * cf - 2.0 * fg * (2.0 * h1 * ksqrtc - h0)) * v_norm;
                    } else if (v_norm > epsilon)
                        integral -= w * M_2_PI * log(v_norm) * cf;
                }
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
    size_t i, j, I, J, N = GaussQR.n;
    double t, s, w, fg, fg_arc, v_norm;
    complex_t cf;
    for (j = 0; j < Qtrial; ++j) {
        for (i = 0; i < Qtrial; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (I = 0; I < N; ++I) {
                for (J = 0; J < N; ++J) {
                    t = GaussQR.x(J), s = GaussQR.x(I), w = GaussQR.w(I) * GaussQR.w(J);
                    fg = trial_space.evaluateShapeFunction(i, s) * trial_space.evaluateShapeFunction(j, t);
                    fg_arc = trial_space.evaluateShapeFunctionDot_01(i, s) * trial_space.evaluateShapeFunctionDot_01(j, t);
                    const complex_t &tangent_p = m_tangent_p(I, J), &tangent = m_tangent(I, J), &h0 = m_h0(I, J), &h1 = m_h1(I, J);
                    fg *= tangent.imag() * tangent_p.imag() + tangent.real() * tangent_p.real();
                    cf = fg_arc - kkc * fg, v_norm = m_v_norm(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += w * h0 * cf;
                        if (der > 0)
                            integral_der += w * (h1 * v_norm * cf + 2.0 * h0 * ksqrtc * fg);
                        if (der > 1)
                            integral_der2 += w * ((h0 * v_norm - h1) * cf - 2.0 * fg * (2.0 * h1 * ksqrtc - h0)) * v_norm;
                    } else if (v_norm > epsilon)
                        integral -= w * M_2_PI * log(v_norm) * cf;
                }
            }
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
    size_t i, j, I, J, N = CGaussQR.n;
    double t, s, w, fg, cf, v_norm;
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < 2 * N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = CGaussQR.x(J % N), s = t * (1. - CGaussQR.x(I)), w = t * CGaussQR.w(I) * CGaussQR.w(J % N);
                    if (J >= N)
                        std::swap(s, t);
                    fg = test_space.evaluateShapeFunction(i, s) * test_space.evaluateShapeFunction(j, t);
                    cf = m_tangent_p_norm(I, J) * m_tangent_norm(I, J) * fg * w, v_norm = m_v_norm(I, J);
                    const complex_t &h0 = m_h0(I, J), &h1 = m_h1(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += h0 * cf;
                        if (der > 0)
                            integral_der += h1 * v_norm * cf;
                        if (der > 1)
                            integral_der2 += (h1 / ksqrtc - h0 * v_norm) * v_norm * cf;
                    } else if (v_norm > epsilon)
                        integral -= M_2_PI * log(v_norm) * cf;
                }
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
    size_t i, j, I, J, N = CGaussQR.n;
    double t, s, w, fg, cf, v_norm;
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < 2 * N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = CGaussQR.x(I), s = t * CGaussQR.x(J % N), w = t * CGaussQR.w(I) * CGaussQR.w(J % N);
                    if (J >= N)
                        std::swap(s, t);
                    fg = swap ? test_space.evaluateShapeFunction(i, s) * test_space.evaluateShapeFunction_01_swapped(j, t)
                              : test_space.evaluateShapeFunction_01_swapped(i, s) * test_space.evaluateShapeFunction(j, t);
                    cf = m_tangent_p_norm(I, J) * m_tangent_norm(I, J) * fg * w, v_norm = m_v_norm(I, J);
                    const complex_t &h0 = m_h0(I, J), &h1 = m_h1(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += h0 * cf;
                        if (der > 0)
                            integral_der += h1 * v_norm * cf;
                        if (der > 1)
                            integral_der2 += (h1 / ksqrtc - h0 * v_norm) * v_norm * cf;
                    } else if (v_norm > epsilon)
                        integral -= M_2_PI * log(v_norm) * cf;
                }
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
    size_t i, j, I, J, N = GaussQR.n;
    double t, s, w, fg, cf, v_norm;
    for (j = 0; j < Qtest; ++j) {
        for (i = 0; i < Qtest; ++i) {
            complex_t integral = czero, integral_der = czero, integral_der2 = czero;
            for (J = 0; J < N; ++J) {
                for (I = 0; I < N; ++I) {
                    t = GaussQR.x(J), s = GaussQR.x(I), w = GaussQR.w(I) * GaussQR.w(J);
                    fg = test_space.evaluateShapeFunction(i, s) * test_space.evaluateShapeFunction(j, t);
                    cf = m_tangent_p_norm(I, J) * m_tangent_norm(I, J) * fg * w, v_norm = m_v_norm(I, J);
                    const complex_t &h0 = m_h0(I, J), &h1 = m_h1(I, J);
                    if (ksqrtca * v_norm > epsilon) {
                        integral += h0 * cf;
                        if (der > 0)
                            integral_der += h1 * v_norm * cf;
                        if (der > 1)
                            integral_der2 += (h1 / ksqrtc - h0 * v_norm) * v_norm * cf;
                    } else if (v_norm > epsilon)
                        integral -= M_2_PI * log(v_norm) * cf;
                }
            }
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
                compute_coinciding(pi, pj);
                double_layer_coinciding(der);
            } else if (is_adjacent(pi, pj, swap)) {
                // adjacent panels
                compute_adjacent(pi, pj, swap);
                double_layer_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(pi, pj);
                double_layer_general(der);
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
        for (j = 0; j < dim_trial; ++j) {
            const auto &pj = *panels[j];
            if (i == j) {
                // coinciding panels
                compute_coinciding(pi, pj);
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
        for (j = 0; j < dim_test; ++j) {
            const auto &pj = *panels[j];
            if (i == j) {
                // coinciding panels
                compute_coinciding(pi, pj);
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
    bool swap;
    // Panel oriented assembly
    for (i = 0; i < dim; ++i) {
        const auto &pi = *panels[i];
        for (j = 0; j < dim; ++j) {
            const auto &pj = *panels[j];
            if (i == j) {
                // coinciding panels
                compute_coinciding(pi, pj);
                double_layer_coinciding(der);
                hypersingular_coinciding(der);
                single_layer_coinciding(der);
            } else if (is_adjacent(pi, pj, swap)) {
                // adjacent panels
                compute_adjacent(pi, pj, swap);
                double_layer_adjacent(swap, der);
                hypersingular_adjacent(swap, der);
                single_layer_adjacent(swap, der);
            } else {
                // disjoint panels
                compute_general(pi, pj);
                double_layer_general(der);
                hypersingular_general(der);
                single_layer_general(der);
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
        }
    }
    double_layer_matrix *= 0.25;
    hypersingular_matrix *= 0.25;
    single_layer_matrix *= 0.25;
    if (der > 1) {
        double_layer_der_matrix *= 0.25;
        hypersingular_der_matrix *= 0.25;
        single_layer_der_matrix *= 0.25;
    }
    if (der > 2) {
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
