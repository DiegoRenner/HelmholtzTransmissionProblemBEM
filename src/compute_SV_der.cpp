#include "mass_matrix.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "single_layer.hpp"
#include "single_layer_der.hpp"
#include "double_layer.hpp"
#include "double_layer_der.hpp"
#include "hypersingular.hpp"
#include "hypersingular_der.hpp"
#include "parametrized_mesh.hpp"
#include "compute_SV_der.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

namespace parametricbem2d {
    typedef std::complex<double> complex_t;
    complex_t ii = complex_t(0, 1);
    complex_t compute_SV_der_simple(const ParametrizedMesh &mesh,
                          unsigned order,
                          const double k_o,
                          const double k_i){
        int numpanels = mesh.getNumPanels();
        // Same trial and test spaces
        DiscontinuousSpace<0> discont_space;
        //DiscontinuousSpace<0> test_space;
        // Space used for interpolation of Dirichlet data
        ContinuousSpace<1> cont_space;
        // Computing V matrix
        std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
        Eigen::MatrixXcd K_o =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_o);
        Eigen::MatrixXcd K_o_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, discont_space, order, k_o);
        std::cout << "double layer helmholtz rhs computed" << std::endl;
        Eigen::MatrixXcd K_i =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_i);
        Eigen::MatrixXcd K_i_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, discont_space, order, k_i);
        std::cout << "double layer helmholtz computed" << std::endl;
        Eigen::MatrixXcd W_i =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
        Eigen::MatrixXcd W_i_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k_i);
        std::cout << "hypersingular helmholtz computed" << std::endl;
        Eigen::MatrixXcd W_o =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k_o);
        Eigen::MatrixXcd W_o_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order,k_o);
        std::cout << "hypersingular helmholtz rhs computed" << std::endl;
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_o);
        Eigen::MatrixXcd V_o_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, discont_space, order, k_o);
        std::cout << "single layer helmholtz rhs computed" << std::endl;
        Eigen::MatrixXcd V_i =
                single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_i);
        Eigen::MatrixXcd V_i_der =
                single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_i);
        std::cout << "single layer helmholtz computed" << std::endl;
        double h = mesh.getPanels()[0]->length()/6.;
        Eigen::VectorXcd H = Eigen::VectorXcd::Ones(numpanels)*h;
        Eigen::MatrixXcd M_cont= (4*H).asDiagonal();
        M_cont.diagonal(-1) = H.segment(0,numpanels-1);
        M_cont.diagonal(+1) = H.segment(0,numpanels-1);
        Eigen::MatrixXcd M_discont = 6*(H).asDiagonal();
        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
        M.block(0, 0, K_o.rows(), K_o.cols()) = M_cont;
        M.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) = M_discont;
        std::cout << "mass matrix computed" << std::endl;

        Eigen::MatrixXcd T(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
        T.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i);
        T.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o-V_i);
        T.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o-W_i;
        T.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o-K_i).transpose();
        //T += M;
        Eigen::MatrixXcd T_der(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
        T_der.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o_der + K_i_der);
        T_der.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o_der-V_i_der);
        T_der.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o_der-W_i_der;
        T_der.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o_der-K_i_der).transpose();
        Eigen::MatrixXcd A_der(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());


        Eigen::MatrixXcd res(2*numpanels,2*numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        res.setIdentity();
        res = llt.transpositionsP()*res;
        res = llt.matrixU() * res;
        res = llt.vectorD().cwiseSqrt().asDiagonal() * res;
        Eigen::MatrixXcd L = res.transpose();
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr_cont(L);
        Eigen::MatrixXcd R = qr_cont.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr_cont.householderQ();
        Eigen::MatrixXcd Trafo = R.inverse()*Q.transpose();
        Eigen::MatrixXcd T_mod = Trafo*T*Trafo.transpose();
        Eigen::MatrixXcd T_mod_der = Trafo*T_der*Trafo.transpose();


        A_der = T_mod_der.conjugate().transpose()*T_mod + T_mod.conjugate().transpose()*T_mod_der;
        //Eigen::BDCSVD<Eigen::MatrixXcd> svd(T_mod);
        //Eigen::VectorXd SVs = svd.singularValues();
        //unsigned num_SVs = SVs.size();

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(T_mod.conjugate().transpose()*T_mod + T_mod + T_mod.conjugate().transpose());
        std::ofstream filename;
        filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_der.dat", std::ios_base::app);
        complex_t val1;
        complex_t val2;

        for ( int i = 0; i < 2*numpanels; i++ ) {
            Eigen::VectorXcd x = sol.eigenvectors().block(0, i, 2 * numpanels, 1);
            val1 = (x).dot(A_der*x);
            val2 = sol.eigenvalues()[i];
            //std::cout << "-----------------"<< std::endl;
            //std::cout << val << std::endl;
            //std::cout << "-----------------"<< std::endl;
            filename << val1.real() << " " << val2.real() << " " ;
        }
        filename << k_o << std::endl;
        filename.close();
        return val1;
    }

    complex_t compute_SV_der(const ParametrizedMesh &mesh,
                             unsigned order,
                             const double k_o,
                             const double k_i){
        int numpanels = mesh.getNumPanels();
        DiscontinuousSpace<0> discont_space;
        ContinuousSpace<1> cont_space;

        Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh, cont_space, cont_space, order);
        //Eigen::MatrixXcd M_discont = mass_matrix::GalerkinMatrix(mesh, discont_space, discont_space, order);
        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        M.block(0,0,numpanels,numpanels) = M_cont;
        M.block(numpanels,numpanels,numpanels,numpanels) = M_cont;
        //Eigen::MatrixXcd M_double = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        //M_double.block(0,0,2*numpanels,2*numpanels) = M;
        //M_double.block(2*numpanels,2*numpanels,2*numpanels,2*numpanels) = M;
        std::cout << "mass matrix computed" << std::endl;

        Eigen::MatrixXcd res(2*numpanels,2*numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        res.setIdentity();
        res = llt.transpositionsP()*res;
        res = llt.matrixU() * res;
        res = llt.vectorD().cwiseSqrt().asDiagonal() * res;
        Eigen::MatrixXcd L = res.transpose();
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr_cont(L);
        Eigen::MatrixXcd R = qr_cont.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr_cont.householderQ();
        Eigen::MatrixXcd Trafo = R.inverse()*Q.transpose();
        std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
        Eigen::MatrixXcd K_o =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k_o);
        Eigen::MatrixXcd K_o_der =
                -double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k_o);
        std::cout << "double layer helmholtz rhs computed" << std::endl;
        Eigen::MatrixXcd K_i =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k_i);
        Eigen::MatrixXcd K_i_der =
                -double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k_i);
        std::cout << "double layer helmholtz computed" << std::endl;
        Eigen::MatrixXcd W_i =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
        Eigen::MatrixXcd W_i_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k_i);
        std::cout << "hypersingular helmholtz computed" << std::endl;
        Eigen::MatrixXcd W_o =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k_o);
        Eigen::MatrixXcd W_o_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order,k_o);
        std::cout << "hypersingular helmholtz rhs computed" << std::endl;
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_o);
        Eigen::MatrixXcd V_o_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k_o);
        std::cout << "single layer helmholtz rhs computed" << std::endl;
        Eigen::MatrixXcd V_i =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
        Eigen::MatrixXcd V_i_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k_i);
        std::cout << "single layer helmholtz computed" << std::endl;

        Eigen::MatrixXcd T = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i);
        T.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o - V_i);
        T.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o -W_i;
        T.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o - K_i).transpose();
        T += M;
        T = Trafo*T*Trafo.transpose();
        Eigen::MatrixXcd T_der = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T_der.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o_der + K_i_der);
        T_der.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o_der - V_i_der);
        T_der.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o_der - W_i_der;
        T_der.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o_der - K_i_der).transpose();
        T_der = Trafo*T_der*Trafo.transpose();

        Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W.block(0,2*numpanels,2*numpanels,2*numpanels) = T;
        W.block(2*numpanels,0,2*numpanels,2*numpanels) = T.transpose().conjugate();
        Eigen::MatrixXcd W_der = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W_der.block(0,2*numpanels,2*numpanels,2*numpanels) = T_der;
        W_der.block(2*numpanels,0,2*numpanels,2*numpanels) = T_der.transpose().conjugate();
        //W = Trafo*W*Trafo.transpose();
        //W_der = Trafo*W_der*Trafo.transpose();
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(W);

        //std::cout << sol.info() << std::endl;
        //std::cout << sol.eigenvectors()*sol.eigenvalues().asDiagonal()*sol.eigenvectors().inverse() - W<< std::endl;
        std::cout << sol.eigenvalues()[numpanels-1]/ sol.eigenvalues()[0] << std::endl;
        //std::cout << sol.matrixV().block(0,2*numpanels-1,1,2*numpanels) << std::endl;
        //std::cout << "********************************************************"<< std::endl;
        //Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol1(M.inverse()*T);
        //std::cout << sol1.eigenvectors().block(0,0,2*numpanels,1) << std::endl;
        //std::cout << sol.matrixV()*sol.singularValues().asDiagonal()*sol.matrixV().conjugate() - W<< std::endl;
        //std::cout << sol.threshold()<< std::endl;
        complex_t normal;
        std::ofstream filename;
        filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_der.dat", std::ios_base::app);
        filename << std::setprecision(25);
        complex_t val1;
        complex_t val2;
        Eigen::VectorXcd x(4*numpanels);
        for ( int i = 0; i < 4*numpanels; i++ ) {
            x = sol.eigenvectors().block(0, i, 4*numpanels, 1);
            //std::cout << x.transpose() << std::endl;
            //std::cout << "***************" << std::endl;
            normal = (x).dot(x);
            //std::cout << normal << std::endl;
            //std::cout << "***************" << std::endl;
            x = x/sqrt(normal);
            //std::cout << x.transpose() << std::endl;
            //std::cout << "***************" << std::endl;
            val1 = (x.dot(W_der*x));//*2.*M_PI/double(numpanels);
            val2 = sol.eigenvalues()[i];
            //std::cout << val1 << " " << val2 << std::endl;
            filename << val1.real() << " " << val2.real() << " " ;
        }
        filename << k_o << std::endl;
        filename.close();
        return val2;
    }

    complex_t compute_SV_der_debug(const ParametrizedMesh &mesh,
                             unsigned order,
                             const double k_o,
                             const double k_i){
        int numpanels = mesh.getNumPanels();
        DiscontinuousSpace<0> discont_space;
        ContinuousSpace<1> cont_space;
        std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
        //Eigen::MatrixXcd W_o =
        //        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_o);
        //Eigen::MatrixXcd W_o_der =
        //        hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k_o);
        //Eigen::MatrixXcd K_o =
        //        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_o);
        //Eigen::MatrixXcd K_o_der =
        //        -double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k_o);
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_o);
        Eigen::MatrixXcd V_o_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, discont_space, order, k_o);
        Eigen::MatrixXcd M = mass_matrix::GalerkinMatrix(mesh, discont_space, discont_space, order);
        Eigen::MatrixXcd M_double = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        M_double.block(0,0,numpanels,numpanels) = M;
        M_double.block(numpanels,numpanels,numpanels,numpanels) = M;
        std::cout << "mass matrix computed" << std::endl;
        Eigen::MatrixXcd res(numpanels,numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        res.setIdentity();
        res = llt.transpositionsP()*res;
        res = llt.matrixU() * res;
        res = llt.vectorD().cwiseSqrt().asDiagonal() * res;
        Eigen::MatrixXcd L = res.transpose();

        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(L);
        Eigen::MatrixXcd R = qr.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr.householderQ();
        Eigen::MatrixXcd Transform = R.inverse()*Q.transpose();
        V_o += M;
        //V_o_der += M;

        V_o = Transform*V_o*Transform.transpose();
        V_o_der = Transform*V_o_der*Transform.transpose();
        //V_o = Transform*V_o*Transform.transpose();
        //V_o_der = Transform*V_o_der*Transform.transpose();
        //K_o = Transform*K_o*Transform.transpose();
        //K_o_der = Transform*K_o_der*Transform.transpose();
        //W_o = Transform*W_o*Transform.transpose();
        //W_o_der = Transform*W_o_der*Transform.transpose();

       Eigen::MatrixXcd W  = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
       W.block(0,numpanels,numpanels,numpanels) = V_o;
       W.block(numpanels,0,numpanels,numpanels) = V_o.transpose().conjugate();
        Eigen::MatrixXcd W_der = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        W_der.block(0,numpanels,numpanels,numpanels) = V_o_der;
        W_der.block(numpanels,0,numpanels,numpanels) = V_o_der.transpose().conjugate();

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(W);
        complex_t normal;
        std::ofstream filename;
        filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_der_debug.dat", std::ios_base::app);
        complex_t val1;
        complex_t val2;
        Eigen::VectorXcd x(2*numpanels);
        complex_t result = sol.eigenvalues()[0];
        for ( int i = 0; i < 2*numpanels; i++ ) {
            x = sol.eigenvectors().block(0, i, 2*numpanels, 1);
            normal = (x).dot(x);
            x = x/sqrt(normal);
            val1 = (x.dot(W_der*x));
            val2 = sol.eigenvalues()[i];
            //std::cout << (K_o*x-val2*x).norm() << std::endl;
            filename << val1.real() << " " << val2.real() << " " ;
        }
        filename << k_o << std::endl;
        filename.close();
        return result;
    }

    double eval_sv(const ParametrizedMesh &mesh,
                                              unsigned order,
                                              const double k_o,
                                              const double k_i){
        int numpanels = mesh.getNumPanels();
        DiscontinuousSpace<0> discont_space;
        ContinuousSpace<1> cont_space;

        Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh, cont_space, cont_space, order);
        //Eigen::MatrixXcd M_discont = mass_matrix::GalerkinMatrix(mesh, discont_space, discont_space, order);
        Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        M.block(0,0,numpanels,numpanels) = M_cont;
        M.block(numpanels,numpanels,numpanels,numpanels) = M_cont;
        //Eigen::MatrixXcd M_double = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        //M_double.block(0,0,2*numpanels,2*numpanels) = M;
        //M_double.block(2*numpanels,2*numpanels,2*numpanels,2*numpanels) = M;

        Eigen::MatrixXcd res(2*numpanels,2*numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        res.setIdentity();
        res = llt.transpositionsP()*res;
        res = llt.matrixU() * res;
        res = llt.vectorD().cwiseSqrt().asDiagonal() * res;
        Eigen::MatrixXcd L = res.transpose();
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr_cont(L);
        Eigen::MatrixXcd R = qr_cont.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr_cont.householderQ();
        Eigen::MatrixXcd Trafo = R.inverse()*Q.transpose();

        Eigen::MatrixXcd K_o =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k_o);
        Eigen::MatrixXcd K_i =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k_i);
        Eigen::MatrixXcd W_i =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
        Eigen::MatrixXcd W_o =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k_o);
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_o);
        Eigen::MatrixXcd V_i =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);

        Eigen::MatrixXcd T = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o +  K_i);
        T.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o - V_i);
        T.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o - W_i;
        T.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o - K_i).transpose();
        T += M;
        T = Trafo*T*Trafo.transpose();

        Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W.block(0,2*numpanels,2*numpanels,2*numpanels) = T;
        W.block(2*numpanels,0,2*numpanels,2*numpanels) = T.conjugate().transpose();
        //W = Trafo*W*Trafo.transpose();
        Eigen::BDCSVD<Eigen::MatrixXcd> sol(T);
        return (sol.singularValues()[2*numpanels-1])*(sol.singularValues()[2*numpanels-1]);
    }

    double eval_sv_debug(const ParametrizedMesh &mesh,
                   unsigned order,
                   const double k_o,
                   const double k_i){
        int numpanels = mesh.getNumPanels();
        DiscontinuousSpace<0> discont_space;
        ContinuousSpace<1> cont_space;
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_o);
        //Eigen::MatrixXcd K_o =
        //        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space,cont_space, order, k_o);
        //Eigen::MatrixXcd W_o =
        //        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_o);
        Eigen::MatrixXcd M = mass_matrix::GalerkinMatrix(mesh, discont_space, discont_space, order);
        Eigen::MatrixXcd res(numpanels,numpanels);
        Eigen::LDLT<Eigen::MatrixXcd> llt(M);
        res.setIdentity();
        res = llt.transpositionsP()*res;
        res = llt.matrixU() * res;
        res = llt.vectorD().cwiseSqrt().asDiagonal() * res;
        Eigen::MatrixXcd L = res.transpose();
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(L);
        Eigen::MatrixXcd R = qr.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXcd Q = qr.householderQ();
        Eigen::MatrixXcd Transform = R.inverse()*Q.transpose();

        V_o += M;
        V_o = Transform*V_o*Transform.transpose();
        //K_o = Transform*K_o*Transform.transpose();
        //W_o = Transform*W_o*Transform.transpose();
        Eigen::BDCSVD<Eigen::MatrixXcd> sol(V_o);

        //std::cout << SVs << std::endl;
        return sol.singularValues()[numpanels-1];
    }


    // Extrapolation based numerical differentation
// with a posteriori error control
// f: handle of a function defined in a neighbourhood of x ∈ R
// x: point at which approximate derivative is desired
// h0: initial distance from x
// rtol: relative target tolerance, atol: absolute tolerance
    double diffex ( std::function<double(double)> f , const double x , const double h0 ,
                    const double rtol , const double atol ) {
        const unsigned nit = 2; // Maximum depth of extrapolation
        Eigen::VectorXd h(nit);
        h[0] = h0 ; // Widths of difference quotients
        Eigen::VectorXd y(nit); // Approximations returned by difference quotients
        y[0] = (f(x+h0) - f(x-h0))/(2*h0); // Widest difference quotients

        // using Aitken-Neville scheme with x = 0, see Code 5.2.3.10
        for ( unsigned i = 1; i < nit ; ++i ){
            // create data points for extrapolation
            h[i] = h[i-1]/2; // Next width half a big
            y[i] = (f(x+h[i]) - f(x-h[i]))/h(i-1) ;
            // Aitken-Neville update
            for (int k = i-1; k >= 0; --k )
                y[k] = y[k+1] - (y[k+1]-y[k] )*h[i]/(h[i]-h[k]) ;
            // termination of extrapolation when desired tolerance is reached
            const double errest = std::abs(y[1]-y[0]) ; // error indicator
            std::cout << y[0] << std::endl;
            if(errest<rtol*std::abs ( y [ 0 ] ) || errest < atol ) //
                break ;
        }
        return y [ 0 ] ; // Return value extrapolated from largest number ofdifference quotients
    }
} // namespace parametricbem2d
