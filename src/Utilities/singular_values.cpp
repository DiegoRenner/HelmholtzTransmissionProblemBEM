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
#include "derivatives.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

namespace parametricbem2d {
    typedef std::complex<double> complex_t;
    complex_t ii = complex_t(0, 1);
    Eigen::MatrixXd sv_1st_der(const ParametrizedMesh &mesh,
                               unsigned order,
                               const complex_t k,
                               const double c_o,
                               const double c_i){
        // get number of panels in mesh and initialize FEM-spaces
        int numpanels = mesh.getNumPanels();
        DiscontinuousSpace<0> discont_space;
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
        Eigen::MatrixXcd K_o =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_o_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_i =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd K_i_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_o =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k, c_o);
        Eigen::MatrixXcd W_o_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_o_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_i =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd V_i_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_i);

        // build solutions operator and it's derivative, project them
        Eigen::MatrixXcd T = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i);
        T.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o - V_i);
        T.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o -W_i;
        T.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o - K_i).transpose();
        T = T + M;
        T = Trafo*T*Trafo.transpose();
        Eigen::MatrixXcd T_der = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T_der.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o_der + K_i_der);
        T_der.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o_der - V_i_der);
        T_der.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o_der - W_i_der;
        T_der.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o_der - K_i_der).transpose();
        T_der = Trafo*T_der*Trafo.transpose();

        // build Wielandt matrix
        Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W.block(0,2*numpanels,2*numpanels,2*numpanels) = T;
        W.block(2*numpanels,0,2*numpanels,2*numpanels) = T.transpose().conjugate();
        Eigen::MatrixXcd W_der = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W_der.block(0,2*numpanels,2*numpanels,2*numpanels) = T_der;
        W_der.block(2*numpanels,0,2*numpanels,2*numpanels) = T_der.transpose().conjugate();

        // get eigenvalues and eigenvectors
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(W);

        // get positive eigenvalue (corresponding to singular value) and compute derivative
        complex_t normal;
        Eigen::MatrixXd res(2*numpanels,2);
        Eigen::VectorXcd x(4*numpanels);
        for ( int i = 0; i < 2*numpanels; i++ ) {
            if ( sol.eigenvalues()[2*i].real() > 0 ){
                res(i,0) = sol.eigenvalues()[2*i].real();
                x = sol.eigenvectors().block(0, 2*i, 4*numpanels, 1);
                normal = (x).dot(x);
                x = x/sqrt(normal);
                res(i,1) = (x.dot(W_der*x)).real();
            } else{
                res(i,0) = sol.eigenvalues()[2*i+1].real();
                x = sol.eigenvectors().block(0, 2*i+1, 4*numpanels, 1);
                normal = (x).dot(x);
                x = x/sqrt(normal);
                res(i,1) = (x.dot(W_der*x)).real();
            }
        }
        return res;
    }

    Eigen::MatrixXd sv_2nd_der(const ParametrizedMesh &mesh,
                               unsigned order,
                               const complex_t k,
                               const double c_o,
                               const double c_i){
        // get number of panels in mesh and initialize FEM-spaces
        int numpanels = mesh.getNumPanels();
        DiscontinuousSpace<0> discont_space;
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
        Eigen::MatrixXcd K_o =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_o_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_o_der2 =
                double_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
        Eigen::MatrixXcd K_i =
                double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd K_i_der =
                double_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd K_i_der2 =
                double_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_i_der2 =
                hypersingular_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd W_o =
                hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k, c_o);
        Eigen::MatrixXcd W_o_der =
                hypersingular_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd W_o_der2 =
                hypersingular_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_o =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_o_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_o_der2 =
                single_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_o);
        Eigen::MatrixXcd V_i =
                single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd V_i_der =
                single_layer_helmholtz_der::GalerkinMatrix(mesh, cont_space, order, k, c_i);
        Eigen::MatrixXcd V_i_der2 =
                single_layer_helmholtz_der2::GalerkinMatrix(mesh, cont_space, order, k, c_i);

        // build solutions operator and it's derivative, project them
        Eigen::MatrixXcd T = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i);
        T.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o - V_i);
        T.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o -W_i;
        T.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o - K_i).transpose();
        T = T + M;
        T = Trafo*T*Trafo.transpose();
        Eigen::MatrixXcd T_der = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T_der.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o_der + K_i_der);
        T_der.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o_der - V_i_der);
        T_der.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o_der - W_i_der;
        T_der.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o_der - K_i_der).transpose();
        T_der = Trafo*T_der*Trafo.transpose();
        Eigen::MatrixXcd T_der2 = Eigen::MatrixXcd::Zero(2*numpanels,2*numpanels);
        T_der2.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o_der2 + K_i_der2);
        T_der2.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o_der2 - V_i_der2);
        T_der2.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o_der2 - W_i_der2;
        T_der2.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                (K_o_der2 - K_i_der2).transpose();
        T_der2 = Trafo*T_der2*Trafo.transpose();

        // build Wielandt matrix
        Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W.block(0,2*numpanels,2*numpanels,2*numpanels) = T;
        W.block(2*numpanels,0,2*numpanels,2*numpanels) = T.transpose().conjugate();
        Eigen::MatrixXcd W_der = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W_der.block(0,2*numpanels,2*numpanels,2*numpanels) = T_der;
        W_der.block(2*numpanels,0,2*numpanels,2*numpanels) = T_der.transpose().conjugate();
        Eigen::MatrixXcd W_der2 = Eigen::MatrixXcd::Zero(4*numpanels,4*numpanels);
        W_der2.block(0,2*numpanels,2*numpanels,2*numpanels) = T_der2;
        W_der2.block(2*numpanels,0,2*numpanels,2*numpanels) = T_der2.transpose().conjugate();
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(W);

        // get positive eigenvalue (corresponding to singular value) and compute derivative
        double temp = 0;
        int m = 5;
        //for(int i = 0; i< 4*numpanels; i++){
        //   if (abs(sol.eigenvectors().col(0)[i])>temp){
        //       temp =abs(sol.eigenvectors().col(0)[i]);
        //       m = i;
        //   }
        //}
        Eigen::MatrixXd res(2*numpanels,3);
        Eigen::MatrixXcd B(4*numpanels,4*numpanels);
        Eigen::VectorXcd r(4*numpanels);
        Eigen::VectorXcd s(4*numpanels);
        Eigen::VectorXcd u_der(4*numpanels);
        Eigen::VectorXcd u_der_temp(4*numpanels);
        Eigen::VectorXcd u_der2(4*numpanels);
        for ( int i = 0; i < 2*numpanels; i++ ) {
            if ( sol.eigenvalues()[2*i].real() > 0 ){
                complex_t rescale = sol.eigenvectors().coeff(m-1,2*i);
                Eigen::MatrixXcd u = sol.eigenvectors().block(0,2*i,4*numpanels,1)/rescale;
                B.block(0,0,4*numpanels,m-1)
                        = (W-sol.eigenvalues()[2*i]*Eigen::MatrixXcd::Identity(4*numpanels,4*numpanels)).block(0,0,4*numpanels,m-1);
                B.block(0,m-1,4*numpanels,4*numpanels-m)
                        = (W-sol.eigenvalues()[2*i]*Eigen::MatrixXcd::Identity(4*numpanels,4*numpanels)).block(0,m,4*numpanels,4*numpanels-m);
                B.block(0,4*numpanels-1,4*numpanels,1) = -u;
                r = -W_der*u;
                u_der_temp = B.inverse()*r;
                complex_t ev_der = u_der_temp[4*numpanels-1];
                u_der.segment(0,m-1) = u_der_temp.segment(0,m-1);
                u_der[m-1] = 0;
                u_der.segment(m,4*numpanels-m) = u_der_temp.segment(m-1,4*numpanels-m);
                s = -W_der2*u -2.*(W_der-ev_der*Eigen::MatrixXcd::Identity(4*numpanels,4*numpanels))*u_der;
                u_der2 = B.inverse()*s;
                complex_t ev_der2 = u_der2[4*numpanels-1];
                res(i,0) = sol.eigenvalues()[2*i].real();
                res(i,1) = ev_der.real();
                res(i,2) = ev_der2.real();
            } else{
                complex_t rescale = sol.eigenvectors().coeff(m-1,2*i+1);
                Eigen::MatrixXcd u = sol.eigenvectors().block(0,2*i+1,4*numpanels,1)/rescale;
                B.block(0,0,4*numpanels,m-1)
                        = (W-sol.eigenvalues()[2*i+1]*Eigen::MatrixXcd::Identity(4*numpanels,4*numpanels)).block(0,0,4*numpanels,m-1);
                B.block(0,m-1,4*numpanels,4*numpanels-m)
                        = (W-sol.eigenvalues()[2*i+1]*Eigen::MatrixXcd::Identity(4*numpanels,4*numpanels)).block(0,m,4*numpanels,4*numpanels-m);
                B.block(0,4*numpanels-1,4*numpanels,1) = -u;
                r = -W_der*u;
                u_der_temp = B.inverse()*r;
                complex_t ev_der = u_der_temp[4*numpanels-1];
                u_der.segment(0,m-1) = u_der_temp.segment(0,m-1);
                u_der[m-1] = 0;
                u_der.segment(m,4*numpanels-m) = u_der_temp.segment(m-1,4*numpanels-m);
                s = -W_der2*u -2.*(W_der-ev_der*Eigen::MatrixXcd::Identity(4*numpanels,4*numpanels))*u_der;
                u_der2 = B.inverse()*s;
                complex_t ev_der2 = u_der2[4*numpanels-1];
                res(i,0) = sol.eigenvalues()[2*i].real();
                res(i,1) = ev_der.real();
                res(i,2) = ev_der2.real();
                res(i,0) = sol.eigenvalues()[2*i+1].real();
            }
        }
        return res;
    }

    double der_by_ext( std::function<double(double)> f ,
                       const double x ,
                       const double h0 ,
                       const double rtol ,
                       const double atol ) {
        const unsigned nit = 2; // Maximum depth of extrapolation
        Eigen::VectorXd h(nit);
        h[0] = h0 ; // Widths of difference quotients
        Eigen::VectorXd y(nit); // Approximations returned by difference quotients
        y[0] = (f(x+h0) - f(x-h0))/(2*h0); // Widest difference quotients

        // using Aitken-Neville scheme with x = 0, see Code 5.2.3.10
        for ( unsigned i = 1; i < nit ; ++i ){
            // create data points for extrapolation
            h[i] = h[i-1]/2; // Next width half as big
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
