#include "singular_values.hpp"
#include <iostream>

typedef std::complex<double> complex_t;
//complex_t ii = complex_t(0, 1);

Eigen::VectorXd sv(const Eigen::MatrixXcd &T,
                   const double *list,
                   const unsigned count){
    // define number of columns for convenience
    unsigned N = T.cols();
    // build Wielandt matrix
    Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(2*N,2*N);
    W.block(0,N,N,N) = T;
    W.block(N,0,N,N) = T.transpose().conjugate();
    // compute eigenvalues of Wielandt matrix
    // which are equivalent to T's singular values
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(W);
    // write out chosen singular values and return
    Eigen::VectorXd res(count);
    unsigned j;
    for (unsigned i = 0; i < count; i++){
        j = list[i];
        if ( sol.eigenvalues()[2*j].real() > 0 ){
            res(i) = sol.eigenvalues()[2*j].real();
        } else{
            res(i) = -sol.eigenvalues()[2*j].real();
        }
    }
    return res;
}

Eigen::MatrixXd sv_1st_der(const Eigen::MatrixXcd &T,
                           const Eigen::MatrixXcd &T_der,
                           const double *list,
                           const unsigned count){
    // get dimensions of operator
    const unsigned N = T.cols();
    // build Wielandt matrix
    Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(2*N,2*N);
    W.block(0,N,N,N) = T;
    W.block(N,0,N,N) = T.transpose().conjugate();
    Eigen::MatrixXcd W_der = Eigen::MatrixXcd::Zero(2*N,2*N);
    W_der.block(0,N,N,N) = T_der;
    W_der.block(N,0,N,N) = T_der.transpose().conjugate();
    // get eigenvalues and eigenvectors
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(W);
    // get positive eigenvalue (corresponding to singular value and compute derivative
    complex_t normal;
    Eigen::MatrixXd res(count,2);
    Eigen::VectorXcd x(2*N);
    unsigned j;
    for ( unsigned i = 0; i < count; i++ ) {
        j = list[i];
        x = sol.eigenvectors().block(0, 2*j, 2*N, 1);
        normal = (x).dot(x);
        x = x/sqrt(normal);
        if ( sol.eigenvalues()[2*j].real() > 0 ){
            res(i,0) = sol.eigenvalues()[2*j].real();
            res(i,1) = (x.dot(W_der*x)).real();
        } else {
            res(i, 0) = -sol.eigenvalues()[2 * j].real();
            res(i, 1) = -(x.dot(W_der * x)).real();
        }
    }
    return res;
}

Eigen::MatrixXd sv_2nd_der(const Eigen::MatrixXcd &T,
                           const Eigen::MatrixXcd &T_der,
                           const Eigen::MatrixXcd &T_der2,
                           const double *list,
                           const unsigned count){
    // get dimensions of operator
    const unsigned N = T.cols();
    // build Wielandt matrix
    Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(2*N,2*N);
    W.block(0,N,N,N) = T;
    W.block(N,0,N,N) = T.transpose().conjugate();
    Eigen::MatrixXcd W_der = Eigen::MatrixXcd::Zero(2*N,2*N);
    W_der.block(0,N,N,N) = T_der;
    W_der.block(N,0,N,N) = T_der.transpose().conjugate();
    Eigen::MatrixXcd W_der2 = Eigen::MatrixXcd::Zero(2*N,2*N);
    W_der2.block(0,N,N,N) = T_der2;
    W_der2.block(N,0,N,N) = T_der2.transpose().conjugate();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> sol(W);
    // get positive eigenvalue (corresponding to singular value) and compute derivative
    Eigen::MatrixXd res(count,3);
    Eigen::MatrixXcd B(2*N,2*N);
    Eigen::VectorXcd r(2*N);
    Eigen::VectorXcd s(2*N);
    Eigen::VectorXcd u_der(2*N);
    Eigen::VectorXcd u_der_temp(2*N);
    Eigen::VectorXcd u_der2(2*N);
    unsigned j;
    for ( unsigned i = 0; i < count; i++ ) {
        // save index of singular value to compute
        j = list[i];

        // choose which entry of eigenvector to normalize
        double temp = 0;
        int m = 5;
        for(unsigned l = 0; l< 2*N; l++){
            if (abs(sol.eigenvectors().col(2*j)[l])>temp){
                temp =abs(sol.eigenvectors().col(0)[l]);
                m = l;
            }
        }
        m += 1;
        complex_t rescale = sol.eigenvectors().coeff(m - 1, 2 * j);
        // normalize eigenvector
        Eigen::MatrixXcd u = sol.eigenvectors().block(0, 2 * j, 2*N, 1) / rescale;
        // build matrix with deleted column from Wielandt matrix and eigenvector
        B.block(0, 0, 2*N, m - 1)
                = (W - sol.eigenvalues()[2 * j] * Eigen::MatrixXcd::Identity(2*N, 2*N))
                .block(0, 0, 2*N, m - 1);
        B.block(0, m - 1, 2*N, 2*N - m)
                = (W - sol.eigenvalues()[2 * j] * Eigen::MatrixXcd::Identity(2*N, 2*N))
                .block(0, m, 2*N, 2*N - m);
        B.block(0, 2*N - 1, 2*N, 1) = -u;
        // compute right hand side
        r = -W_der * u;
        // solve linear system of equations for derivative of eigenvalue and eigenvector
        u_der_temp = B.inverse() * r;
        complex_t ev_der = u_der_temp[2*N - 1];
        // construct eigenvector using derivative of normalization condition
        u_der.segment(0, m - 1) = u_der_temp.segment(0, m - 1);
        u_der[m - 1] = 0;
        u_der.segment(m, 2*N - m) = u_der_temp.segment(m - 1, 2*N - m);
        // compute right hand side for second derivative
        s = -W_der2 * u -
            2. * (W_der - ev_der * Eigen::MatrixXcd::Identity(2*N, 2*N)) * u_der;
        // solve linear system of equations second derivative of eigenvector and eigenvalue
        u_der2 = B.inverse() * s;
        complex_t ev_der2 = u_der2[2*N - 1];
        // check sign of eigenvalue of Wielandt matrix to retrieve correct eigenvalue of T
        if (sol.eigenvalues()[2 * j].real() > 0) {
            res(i, 0) = sol.eigenvalues()[2 * j].real();
            res(i, 1) = ev_der.real();
            res(i, 2) = ev_der2.real();
        } else {
            res(i, 0) = -sol.eigenvalues()[2 * j].real();
            res(i, 1) = -ev_der.real();
            res(i, 2) = -ev_der2.real();
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
        if(errest<rtol*std::abs ( y [ 0 ] ) || errest < atol ) //
            break ;
    }
    return y [ 0 ] ; // Return value extrapolated from largest number ofdifference quotients
}
