#include "gen_sol_op.hpp"
#include "mass_matrix.hpp"
// #include "double_layer.hpp"
// #include "double_layer_der.hpp"
// #include "double_layer_der2.hpp"
// #include "hypersingular.hpp"
// #include "hypersingular_der.hpp"
// #include "hypersingular_der2.hpp"
// #include "single_layer.hpp"
// #include "single_layer_der.hpp"
// #include "single_layer_der2.hpp"
// #include <iostream>

typedef std::complex<double> complex_t;

SolutionsOperator::SolutionsOperator(const ParametrizedMesh &mesh_in,
                                     unsigned int order,
                                     const AbstractBEMSpace &test_space_in,
                                     const AbstractBEMSpace &trial_space_in)
: mesh(mesh_in), test_space(test_space_in), trial_space(trial_space_in)
{
    GaussQR = getGaussQR(order, 0., 1.);
    CGaussQR = getCGaussQR(order);
    size_t numpanels = mesh.getNumPanels();
    dim_test = test_space.getSpaceDim(numpanels);
    dim_trial = trial_space.getSpaceDim(numpanels);
    size_t dim = dim_test + dim_trial;
    // assemble mass matrix
    M = mass_matrix::GalerkinMatrix(mesh, trial_space, test_space, GaussQR);
    // compute matrix for projection onto ortogonal FEM-spaces
    Eigen::MatrixXd A;
    A.setZero(dim, dim);
    A.block(0, 0, dim_test, dim_trial) = M;
    A.block(dim_test, dim_trial, dim_trial, dim_test) = M.transpose().eval();
    Eigen::LDLT<Eigen::MatrixXd> llt(A);
    lu = Eigen::PartialPivLU<Eigen::MatrixXcd>
        (llt.transpositionsP().transpose() * llt.matrixL().toDenseMatrix() * llt.vectorD().cwiseSqrt().asDiagonal());
}

Eigen::MatrixXcd SolutionsOperator::project(const Eigen::MatrixXcd &T) const {
    return lu.solve(lu.solve(T).transpose().eval()).transpose().eval();
}


void SolutionsOperator::gen_sol_op(const complex_t &k, double c_o, double c_i,
                                   Eigen::MatrixXcd &T) const {
    GalerkinMatrixBuilder builder(mesh, test_space, trial_space, GaussQR, CGaussQR);
    gen_sol_op_in(builder, k, c_o, c_i, T);
}

void SolutionsOperator::gen_sol_op(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                                   Eigen::MatrixXcd &T) const {
    gen_sol_op_in(builder, k, c_o, c_i, T);
}

void SolutionsOperator::gen_sol_op_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                                      Eigen::MatrixXcd &T) const {
    T.resize(dim_test + dim_trial, dim_trial + dim_test);
    Eigen::MatrixXcd K_i, K_o, V_i, V_o, W_i, W_o;
    if (builder.testTrialSpacesAreEqual()) {
        builder.assembleAll(k, c_i);
        K_i = builder.getDoubleLayer();
        W_i = builder.getHypersingular();
        V_i = builder.getSingleLayer();
        builder.assembleAll(k, c_o);
        K_o = builder.getDoubleLayer();
        W_o = builder.getHypersingular();
        V_o = builder.getSingleLayer();
    } else {
        builder.assembleDoubleLayer(k, c_i);
        K_i = builder.getDoubleLayer();
        builder.assembleHypersingular(k, c_i);
        W_i = builder.getHypersingular();
        builder.assembleSingleLayer(k, c_i);
        V_i = builder.getSingleLayer();
        builder.assembleDoubleLayer(k, c_o);
        K_o = builder.getDoubleLayer();
        builder.assembleHypersingular(k, c_o);
        W_o = builder.getHypersingular();
        builder.assembleSingleLayer(k, c_o);
        V_o = builder.getSingleLayer();
    }
    // assemble solutions operator matrix
    T.block(0, 0, dim_test, dim_trial) = M - K_o + K_i;
    T.block(dim_test, 0, dim_trial, dim_trial) = W_o - W_i;
    T.block(0, dim_trial, dim_test, dim_test) = V_o - V_i;
    T.block(dim_test, dim_trial, dim_trial, dim_test) = (M + K_o - K_i).transpose();
    // project onto orthogonal FEM spaces
    T = project(T);
}

void SolutionsOperator::gen_sol_op_1st_der(const complex_t &k, double c_o, double c_i,
                                           Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der) const {
    GalerkinMatrixBuilder builder(mesh, test_space, trial_space, GaussQR, CGaussQR);
    gen_sol_op_1st_der_in(builder, k, c_o, c_i, T, T_der);
}

void SolutionsOperator::gen_sol_op_1st_der(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                                           Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der) const {
    gen_sol_op_1st_der_in(builder, k, c_o, c_i, T, T_der);
}

void SolutionsOperator::gen_sol_op_1st_der_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                                              Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der) const {
    T.resize(dim_test + dim_trial, dim_trial + dim_test);
    T_der.resize(dim_test + dim_trial, dim_trial + dim_test);
    Eigen::MatrixXcd K_i, K_o, V_i, V_o, W_i, W_o;
    Eigen::MatrixXcd K_der_i, K_der_o, V_der_i, V_der_o, W_der_i, W_der_o;
    if (builder.testTrialSpacesAreEqual()) {
        builder.assembleAll(k, c_i, 1);
        K_i = builder.getDoubleLayer(0);
        W_i = builder.getHypersingular(0);
        V_i = builder.getSingleLayer(0);
        K_der_i = builder.getDoubleLayer(1);
        W_der_i = builder.getHypersingular(1);
        V_der_i = builder.getSingleLayer(1);
        builder.assembleAll(k, c_o, 1);
        K_o = builder.getDoubleLayer(0);
        W_o = builder.getHypersingular(0);
        V_o = builder.getSingleLayer(0);
        K_der_o = builder.getDoubleLayer(1);
        W_der_o = builder.getHypersingular(1);
        V_der_o = builder.getSingleLayer(1);
    } else {
        builder.assembleDoubleLayer(k, c_i, 1);
        K_i = builder.getDoubleLayer(0);
        K_der_i = builder.getDoubleLayer(1);
        builder.assembleHypersingular(k, c_i, 1);
        W_i = builder.getHypersingular(0);
        W_der_i = builder.getHypersingular(1);
        builder.assembleSingleLayer(k, c_i, 1);
        V_i = builder.getSingleLayer(0);
        V_der_i = builder.getSingleLayer(1);
        builder.assembleDoubleLayer(k, c_o, 1);
        K_o = builder.getDoubleLayer(0);
        K_der_o = builder.getDoubleLayer(1);
        builder.assembleHypersingular(k, c_o, 1);
        W_o = builder.getHypersingular(0);
        W_der_o = builder.getHypersingular(1);
        builder.assembleSingleLayer(k, c_o, 1);
        V_o = builder.getSingleLayer(0);
        V_der_o = builder.getSingleLayer(1);
    }
    // assemble solutions operator matrix and its derivative
    T.block(0, 0, dim_test, dim_trial) = M - K_o + K_i;
    T.block(dim_test, 0, dim_trial, dim_trial) = W_o - W_i;
    T.block(0, dim_trial, dim_test, dim_test) = V_o - V_i;
    T.block(dim_test, dim_trial, dim_trial, dim_test) = (M + K_o - K_i).transpose();
    T_der.block(0, 0, dim_test, dim_trial) = -K_der_o + K_der_i;
    T_der.block(dim_test, 0, dim_trial, dim_trial) = W_der_o - W_der_i;
    T_der.block(0, dim_trial, dim_test, dim_test) = V_der_o - V_der_i;
    T_der.block(dim_test, dim_trial, dim_trial, dim_test) = (K_der_o - K_der_i).transpose();
    // project onto orthogonal FEM spaces
    T = project(T);
    T_der = project(T_der);
}

void SolutionsOperator::gen_sol_op_2nd_der(const complex_t &k, double c_o, double c_i,
                                           Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2) const {
    GalerkinMatrixBuilder builder(mesh, test_space, trial_space, GaussQR, CGaussQR);
    gen_sol_op_2nd_der_in(builder, k, c_o, c_i, T, T_der, T_der2);
}

void SolutionsOperator::gen_sol_op_2nd_der(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                                           Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2) const {
    gen_sol_op_2nd_der_in(builder, k, c_o, c_i, T, T_der, T_der2);
}

void SolutionsOperator::gen_sol_op_2nd_der_in(GalerkinMatrixBuilder &builder, const complex_t &k, double c_o, double c_i,
                                              Eigen::MatrixXcd &T, Eigen::MatrixXcd &T_der, Eigen::MatrixXcd &T_der2) const {
    T.resize(dim_test + dim_trial, dim_trial + dim_test);
    T_der.resize(dim_test + dim_trial, dim_trial + dim_test);
    T_der2.resize(dim_test + dim_trial, dim_trial + dim_test);
    Eigen::MatrixXcd K_i, K_o, V_i, V_o, W_i, W_o;
    Eigen::MatrixXcd K_der_i, K_der_o, V_der_i, V_der_o, W_der_i, W_der_o;
    Eigen::MatrixXcd K_der2_i, K_der2_o, V_der2_i, V_der2_o, W_der2_i, W_der2_o;
    if (builder.testTrialSpacesAreEqual()) {
        builder.assembleAll(k, c_i, 2);
        K_i = builder.getDoubleLayer(0);
        W_i = builder.getHypersingular(0);
        V_i = builder.getSingleLayer(0);
        K_der_i = builder.getDoubleLayer(1);
        W_der_i = builder.getHypersingular(1);
        V_der_i = builder.getSingleLayer(1);
        K_der2_i = builder.getDoubleLayer(2);
        W_der2_i = builder.getHypersingular(2);
        V_der2_i = builder.getSingleLayer(2);
        builder.assembleAll(k, c_o, 2);
        K_o = builder.getDoubleLayer(0);
        W_o = builder.getHypersingular(0);
        V_o = builder.getSingleLayer(0);
        K_der_o = builder.getDoubleLayer(1);
        W_der_o = builder.getHypersingular(1);
        V_der_o = builder.getSingleLayer(1);
        K_der2_o = builder.getDoubleLayer(2);
        W_der2_o = builder.getHypersingular(2);
        V_der2_o = builder.getSingleLayer(2);
    } else {
        builder.assembleDoubleLayer(k, c_i, 2);
        K_i = builder.getDoubleLayer(0);
        K_der_i = builder.getDoubleLayer(1);
        K_der2_i = builder.getDoubleLayer(2);
        builder.assembleHypersingular(k, c_i, 2);
        W_i = builder.getHypersingular(0);
        W_der_i = builder.getHypersingular(1);
        W_der2_i = builder.getHypersingular(2);
        builder.assembleSingleLayer(k, c_i, 2);
        V_i = builder.getSingleLayer(0);
        V_der_i = builder.getSingleLayer(1);
        V_der2_i = builder.getSingleLayer(2);
        builder.assembleDoubleLayer(k, c_o, 2);
        K_o = builder.getDoubleLayer(0);
        K_der_o = builder.getDoubleLayer(1);
        K_der2_o = builder.getDoubleLayer(2);
        builder.assembleHypersingular(k, c_o, 2);
        W_o = builder.getHypersingular(0);
        W_der_o = builder.getHypersingular(1);
        W_der2_o = builder.getHypersingular(2);
        builder.assembleSingleLayer(k, c_o, 2);
        V_o = builder.getSingleLayer(0);
        V_der_o = builder.getSingleLayer(1);
        V_der2_o = builder.getSingleLayer(2);
    }
    // assemble solutions operator matrix and its derivatives
    T.block(0, 0, dim_test, dim_trial) = M - K_o + K_i;
    T.block(dim_test, 0, dim_trial, dim_trial) = W_o - W_i;
    T.block(0, dim_trial, dim_test, dim_test) = V_o - V_i;
    T.block(dim_test, dim_trial, dim_trial, dim_test) = (M + K_o - K_i).transpose();
    T_der.block(0, 0, dim_test, dim_trial) = -K_der_o + K_der_i;
    T_der.block(dim_test, 0, dim_trial, dim_trial) = W_der_o - W_der_i;
    T_der.block(0, dim_trial, dim_test, dim_test) = V_der_o - V_der_i;
    T_der.block(dim_test, dim_trial, dim_trial, dim_test) = (K_der_o - K_der_i).transpose();
    T_der2.block(0, 0, dim_test, dim_trial) = -K_der2_o + K_der2_i;
    T_der2.block(dim_test, 0, dim_trial, dim_trial) = W_der2_o - W_der2_i;
    T_der2.block(0, dim_trial, dim_test, dim_test) = V_der2_o - V_der2_i;
    T_der2.block(dim_test, dim_trial, dim_trial, dim_test) = (K_der2_o - K_der2_i).transpose();
    // project onto orthogonal FEM spaces
    T = project(T);
    T_der = project(T_der);
    T_der2 = project(T_der2);
}
