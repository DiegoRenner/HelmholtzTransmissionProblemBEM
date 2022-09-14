#include "gen_sol_op.hpp"
#include "mass_matrix.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "single_layer.hpp"
#include "single_layer_der.hpp"
#include "single_layer_der2.hpp"
#include "double_layer.hpp"
#include "double_layer_der.hpp"
#include "double_layer_der2.hpp"
#include "hypersingular.hpp"
#include "hypersingular_der.hpp"
#include "hypersingular_der2.hpp"
#include "parametrized_mesh.hpp"


typedef std::complex<double> complex_t;
    Eigen::MatrixXcd gen_sol_op(const ParametrizedMesh &mesh,
                                unsigned order,
                                const complex_t k,
                                const double c_o,
                                const double c_i){
        // get number of panels in mesh and initialize FEM-spaces
        int numpanels = mesh.getNumPanels();
        ContinuousSpace<1> cont_space;

        // compute mass matrix
        Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh, cont_space, cont_space, order);
        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        M.block(0,0,numpanels,numpanels) = M_cont;
        M.block(numpanels,numpanels,numpanels,numpanels) = M_cont;

        // compute matrix for projection onto ortogonal FEM-spaces
        Eigen::MatrixXcd lt(2*numpanels,2*numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        lt.setIdentity();
        lt = llt.transpositionsP()*lt;
        lt = llt.matrixU() * lt;
        lt = llt.vectorD().cwiseSqrt().asDiagonal() * lt;
        Eigen::MatrixXcd L = lt.transpose();
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr_cont(L);
        Eigen::MatrixXcd R = qr_cont.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr_cont.householderQ();
        Eigen::MatrixXcd Trafo = R.inverse()*Q.transpose();

        // compute operator matrices and their derivatives on inner and outer domain
        Eigen::MatrixXcd K_o =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_i =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_o =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k, c_o);
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_i =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);

        // build solutions operator and it's derivative, project them
        Eigen::MatrixXcd T = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T.block(0, 0, numpanels, numpanels) = (-K_o + K_i);
        T.block(0, numpanels, numpanels, numpanels) = (V_o - V_i);
        T.block(numpanels, 0, numpanels, numpanels) = W_o -W_i;
        T.block(numpanels, numpanels, numpanels, numpanels) =
                (K_o - K_i).transpose();
        T = T + M;
        T = Trafo*T*Trafo.transpose();

        return T;
    }

    Eigen::MatrixXcd gen_sol_op_1st_der(const ParametrizedMesh &mesh,
                                unsigned order,
                                const complex_t k,
                                const double c_o,
                                const double c_i){
        // get number of panels in mesh and initialize FEM-spaces
        int numpanels = mesh.getNumPanels();
        ContinuousSpace<1> cont_space;

        // compute mass matrix
        Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh, cont_space, cont_space, order);
        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        M.block(0,0,numpanels,numpanels) = M_cont;
        M.block(numpanels,numpanels,numpanels,numpanels) = M_cont;

        // compute matrix for projection onto ortogonal FEM-sapces
        Eigen::MatrixXcd lt(2*numpanels,2*numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        lt.setIdentity();
        lt = llt.transpositionsP()*lt;
        lt = llt.matrixU() * lt;
        lt = llt.vectorD().cwiseSqrt().asDiagonal() * lt;
        Eigen::MatrixXcd L = lt.transpose();
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr_cont(L);
        Eigen::MatrixXcd R = qr_cont.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr_cont.householderQ();
        Eigen::MatrixXcd Trafo = R.inverse()*Q.transpose();

        // compute operator matrices and their derivatives on inner and outer domain
        Eigen::MatrixXcd K_o_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_i_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_o_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_o_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_i_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_i);

        // build solutions operator and it's derivative, project them
        Eigen::MatrixXcd T_der = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T_der.block(0, 0, numpanels, numpanels) = (-K_o_der + K_i_der);
        T_der.block(0, numpanels, numpanels, numpanels) = (V_o_der - V_i_der);
        T_der.block(numpanels, 0, numpanels, numpanels) = W_o_der - W_i_der;
        T_der.block(numpanels, numpanels, numpanels, numpanels) =
                (K_o_der - K_i_der).transpose();
        T_der = Trafo*T_der*Trafo.transpose();

        return T_der;
    }

    Eigen::MatrixXcd gen_sol_op_2nd_der(const ParametrizedMesh &mesh,
                                        unsigned order,
                                        const complex_t k,
                                        const double c_o,
                                        const double c_i){
        // get number of panels in mesh and initialize FEM-spaces
        int numpanels = mesh.getNumPanels();
        ContinuousSpace<1> cont_space;

        // compute mass matrix
        Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh, cont_space, cont_space, order);
        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        M.block(0,0,numpanels,numpanels) = M_cont;
        M.block(numpanels,numpanels,numpanels,numpanels) = M_cont;

        // compute matrix for projection onto ortogonal FEM-sapces
        Eigen::MatrixXcd lt(2*numpanels,2*numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        lt.setIdentity();
        lt = llt.transpositionsP()*lt;
        lt = llt.matrixU() * lt;
        lt = llt.vectorD().cwiseSqrt().asDiagonal() * lt;
        Eigen::MatrixXcd L = lt.transpose();
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr_cont(L);
        Eigen::MatrixXcd R = qr_cont.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr_cont.householderQ();
        Eigen::MatrixXcd Trafo = R.inverse()*Q.transpose();

        // compute operator matrices and their derivatives on inner and outer domain
        Eigen::MatrixXcd K_o_der2 =
                double_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_i_der2 =
                double_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i_der2 =
                hypersingular_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_o_der2 =
                hypersingular_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_o_der2 =
                single_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_i_der2 =
                single_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_i);

        // build solutions operator and it's derivative, project them
        Eigen::MatrixXcd T_der2 = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T_der2.block(0, 0, numpanels, numpanels) = (-K_o_der2 + K_i_der2);
        T_der2.block(0, numpanels, numpanels, numpanels) = (V_o_der2 - V_i_der2);
        T_der2.block(numpanels, 0, numpanels, numpanels) = W_o_der2 - W_i_der2;
        T_der2.block(numpanels, numpanels, numpanels, numpanels) =
                (K_o_der2 - K_i_der2).transpose();
        T_der2 = Trafo*T_der2*Trafo.transpose();

        return T_der2;
   }
