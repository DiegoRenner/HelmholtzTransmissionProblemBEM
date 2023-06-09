#include "gen_sol_op.hpp"
#include "mass_matrix.hpp"
#include "galerkin_all.hpp"
#include "galerkin_all_der.hpp"
#include "single_layer.hpp"
#include "single_layer_der.hpp"
#include "single_layer_der2.hpp"
#include "double_layer.hpp"
#include "double_layer_der.hpp"
#include "double_layer_der2.hpp"
#include "hypersingular.hpp"
#include "hypersingular_der.hpp"
#include "hypersingular_der2.hpp"

typedef std::complex<double> complex_t;

void SolutionsOperator::tridiagonal_ldl(const Eigen::MatrixXcd &A, Eigen::VectorXcd &l, Eigen::VectorXcd &d) const {
    d(0) = A(0, 0);
    unsigned int n = A.rows(), k = 1;
    for (; k < n; ++k) {
        const complex_t &a = A(k, k - 1);
        d(k) = A(k, k) - (l(k - 1) = a / d(k - 1)) * a;
    }
}

Eigen::MatrixXcd SolutionsOperator::tridiagonal_lu(const Eigen::MatrixXcd &B) const {
    unsigned int n = B.rows(), k = 1;
    Eigen::MatrixXcd res(n, n);
    res.row(0) = B.row(0);
    for (; k < n; ++k) {
        res.row(k) = B.row(k) - L(k - 1) * res.row(k - 1);
        res.row(k - 1) /= U(k - 1);
    }
    res.row(n - 1) /= U(n - 1);
    return res;
}

SolutionsOperator::SolutionsOperator(const ParametrizedMesh &mesh_in, unsigned order_in) : mesh(mesh_in) {
    order = order_in;
    GaussQR = getGaussQR(order, 0., 1.);
    CGaussQR = getCGaussQR(order);
    numpanels = mesh.getNumPanels();
    L.resize(2 * numpanels - 1);
    U.resize(2 * numpanels);
    // compute mass matrix
    M = mass_matrix::GalerkinMatrix(mesh, cont_space, cont_space, GaussQR);
    // compute matrix for projection onto ortogonal FEM-spaces
    Eigen::VectorXcd l(numpanels - 1), u(numpanels);
    tridiagonal_ldl(M, l, u);
    u = u.cwiseSqrt();
    L << l, 0.0, l;
    U << u, u;
}

Eigen::MatrixXcd SolutionsOperator::gen_sol_op(const complex_t &k, double c_o, double c_i) const {
    Eigen::MatrixXcd T;
    T.setZero(2 * numpanels, 2 * numpanels);
    auto ul = T.block(0, 0, numpanels, numpanels);
    auto ur = T.block(0, numpanels, numpanels, numpanels);
    auto dl = T.block(numpanels, 0, numpanels, numpanels);
    auto dr = T.block(numpanels, numpanels, numpanels, numpanels);
    GalerkinMatrices(mesh, cont_space, GaussQR, CGaussQR, k, c_i, c_o, ul, ur, dl);
    dr = M + ul.transpose();
    ul = M - ul;
    return tridiagonal_lu(tridiagonal_lu(T).transpose()).transpose();
}

std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> SolutionsOperator::gen_sol_op_with_1st_der(const complex_t &k, double c_o, double c_i) const {
    Eigen::MatrixXcd T, T_der;
    T.setZero(2 * numpanels, 2 * numpanels);
    T_der.setZero(2 * numpanels, 2 * numpanels);
    auto ul = T.block(0, 0, numpanels, numpanels);
    auto ur = T.block(0, numpanels, numpanels, numpanels);
    auto dl = T.block(numpanels, 0, numpanels, numpanels);
    auto dr = T.block(numpanels, numpanels, numpanels, numpanels);
    auto ul_der = T_der.block(0, 0, numpanels, numpanels);
    auto ur_der = T_der.block(0, numpanels, numpanels, numpanels);
    auto dl_der = T_der.block(numpanels, 0, numpanels, numpanels);
    auto dr_der = T_der.block(numpanels, numpanels, numpanels, numpanels);
    GalerkinMatrices_der(mesh, cont_space, GaussQR, CGaussQR, k, c_i, c_o, ul, ul_der, ur, ur_der, dl, dl_der);
    dr = M + ul.transpose();
    ul = M - ul;
    dr_der = ul_der.transpose();
    ul_der = -ul_der;
    return std::make_pair(tridiagonal_lu(tridiagonal_lu(T).transpose()).transpose(), tridiagonal_lu(tridiagonal_lu(T_der).transpose()).transpose());
}

Eigen::MatrixXcd SolutionsOperator::gen_sol_op_1st_der(const complex_t &k, double c_o, double c_i) const {
    // compute operator matrices and their derivatives on inner and outer domain
    Eigen::MatrixXcd K_der = double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, GaussQR, CGaussQR, k, c_i, c_o);
    Eigen::MatrixXcd W_der = hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, GaussQR, CGaussQR, k, c_i, c_o);
    Eigen::MatrixXcd V_der = single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, GaussQR, CGaussQR, k, c_i, c_o);
    // build solutions operator and it's derivative, project them
    Eigen::MatrixXcd T_der = Eigen::MatrixXcd::Zero(2 * numpanels, 2 * numpanels);
    T_der.block(0, 0, numpanels, numpanels) = -K_der;
    T_der.block(0, numpanels, numpanels, numpanels) = V_der;
    T_der.block(numpanels, 0, numpanels, numpanels) = W_der;
    T_der.block(numpanels, numpanels, numpanels, numpanels) = K_der.transpose();
    return tridiagonal_lu(tridiagonal_lu(T_der).transpose()).transpose();
}

Eigen::MatrixXcd SolutionsOperator::gen_sol_op_2nd_der(const complex_t &k, double c_o, double c_i) const {
    // compute operator matrices and their derivatives on inner and outer domain
    Eigen::MatrixXcd K_der2 = double_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, cont_space, GaussQR, CGaussQR, k, c_i, c_o);
    Eigen::MatrixXcd W_der2 = hypersingular_helmholtz_der2::GalerkinMatrix(mesh, cont_space, GaussQR, CGaussQR, k, c_i, c_o);
    Eigen::MatrixXcd V_der2 = single_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, GaussQR, CGaussQR, k, c_i, c_o);
    // build solutions operator and it's derivative, project them
    Eigen::MatrixXcd T_der2 = Eigen::MatrixXcd::Zero(2 * numpanels, 2 * numpanels);
    T_der2.block(0, 0, numpanels, numpanels) = -K_der2;
    T_der2.block(0, numpanels, numpanels, numpanels) = V_der2;
    T_der2.block(numpanels, 0, numpanels, numpanels) = W_der2;
    T_der2.block(numpanels, numpanels, numpanels, numpanels) = K_der2.transpose();
    return tridiagonal_lu(tridiagonal_lu(T_der2).transpose()).transpose();
}

Eigen::MatrixXcd SolutionsOperator::K_cont(const std::complex<double> &k, double c) const {
    return double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, GaussQR, CGaussQR, k, 0., c);
}

Eigen::MatrixXcd SolutionsOperator::K_cont_discont(const std::complex<double> &k, double c) const {
    return double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, GaussQR, CGaussQR, k, 0., c);
}

Eigen::MatrixXcd SolutionsOperator::K_discont_cont(const std::complex<double> &k, double c) const {
    return double_layer_helmholtz::GalerkinMatrix(mesh, discont_space, cont_space, GaussQR, CGaussQR, k, 0., c);
}

Eigen::MatrixXcd SolutionsOperator::V_cont(const std::complex<double> &k, double c) const {
    return single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, GaussQR, CGaussQR, k, 0., c);
}

Eigen::MatrixXcd SolutionsOperator::V_discont(const std::complex<double> &k, double c) const {
    return single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, GaussQR, CGaussQR, k, 0., c);
}

Eigen::MatrixXcd SolutionsOperator::W_cont(const std::complex<double> &k, double c) const {
    return hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, GaussQR, CGaussQR, k, 0., c);
}
