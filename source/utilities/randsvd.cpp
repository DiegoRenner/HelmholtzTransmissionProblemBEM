/**
 * \file randsvd.cpp
 *
 * \brief This file contains implementation of routines for creating
 * random Gaussian matrices and for approximating the smallest
 * singular value of a given invertible matrix.
 *
 * This file is a part of the HelmholtzTransmissionProblemBEM library.
 *
 * (c) 2023 Luka Marohnić
 */

#include "randsvd.hpp"
#include <complex>
#include <random>
#include <iostream>

namespace randomized_svd {

    Eigen::MatrixXcd randGaussian(int nr, int nc) {
        std::random_device rd {};
        std::mt19937 gen { rd() };
        std::uniform_real_distribution<> d { 0, 1 };
        Eigen::MatrixXcd R = Eigen::MatrixXcd::Zero(nr, nc);
        for (int i = 0; i < nr; ++i) for (int j = 0; j < nc; ++j) {
            R(i, j) = std::complex<double>(d(gen), d(gen));
        }
        return R * M_SQRT1_2;
    }

    double sv(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &R, int q) {
        int nr = R.rows(), nc = R.cols();
        Eigen::MatrixXcd Q, thinQ;
        thinQ.setIdentity(nr, nc);
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu_decomp(T);
        Q = lu_decomp.solve(R).colPivHouseholderQr().matrixQ() * thinQ;
        for (int i = 0; i < q; ++i) {
            Q = lu_decomp.adjoint().solve(Q).colPivHouseholderQr().matrixQ() * thinQ;
            Q = lu_decomp.solve(Q).colPivHouseholderQr().matrixQ() * thinQ;
        }
        auto svd = lu_decomp.adjoint().solve(Q).bdcSvd();
        return 1.0 / svd.singularValues()(0);
    }

    Eigen::Vector2d sv_der(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &T_der, const Eigen::MatrixXcd &R, int q) {
        unsigned int N = R.rows();
        Eigen::Vector2d res;
        Eigen::MatrixXcd W, W_der;
        W.setZero(2 * N, 2 * N);
        W_der.setZero(2 * N, 2 * N);
        W.block(0, N, N, N) = T;
        W.block(N, 0, N, N) = T.adjoint();
        W_der.block(0, N, N, N) = T_der;
        W_der.block(N, 0, N, N) = T_der.adjoint();
        // the smallest singular value
        double s = sv(T, R, q);
        // get the corresponding eigenvector of the Wielandt
        // matrix W as the last column in the matrix Q of
        // a QR factorization of W - s * I
        Eigen::MatrixXcd A = W - (s * Eigen::VectorXcd::Ones(2 * N)).asDiagonal().toDenseMatrix();
        Eigen::MatrixXcd Q = A.colPivHouseholderQr().matrixQ();
        Eigen::VectorXcd x = Q.col(2 * N - 1);
        // compute the derivative of s
        res(0) = s;
        res(1) = x.dot(W_der * x).real();
        return res;
    }

    Eigen::Vector3d sv_der2(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &T_der, const Eigen::MatrixXcd T_der2, const Eigen::MatrixXcd &R, int q) {
        unsigned int N = R.rows();
        Eigen::Vector3d res;
        Eigen::MatrixXcd W, W_der, W_der2, B(2 * N, 2 * N);
        W.setZero(2 * N, 2 * N);
        W_der.setZero(2 * N, 2 * N);
        W_der2.setZero(2 * N, 2 * N);
        W.block(0, N, N, N) = T;
        W.block(N, 0, N, N) = T.adjoint();
        W_der.block(0, N, N, N) = T_der;
        W_der.block(N, 0, N, N) = T_der.adjoint();
        W_der2.block(0, N, N, N) = T_der2;
        W_der2.block(N, 0, N, N) = T_der2.adjoint();
        // smallest singular value
        double s = sv(T, R, q);
        // compute the first derivative of s
        Eigen::MatrixXcd A = W - (s * Eigen::VectorXcd::Ones(2 * N)).asDiagonal().toDenseMatrix();
        auto qr = A.colPivHouseholderQr();
        Eigen::MatrixXcd Q = qr.matrixQ();
        Eigen::VectorXcd x = Q.col(2 * N - 1), u = x;
        x.normalize();
        res(0) = s;
        res(1) = x.dot(W_der * x).real();
        // compute the second derivative of s
        double temp = 0;
        int m = 5;
        for (unsigned l = 0; l < 2 * N; l++) {
            if (abs(u.coeff(l)) > temp) {
                temp = abs(u.coeff(l));
                m = l;
            }
        }
        m += 1;
        // normalize eigenvector
        u /= u.coeff(m - 1);
        // build matrix with deleted column from Wielandt matrix and eigenvector
        W -= (s * Eigen::VectorXcd::Ones(2 * N)).asDiagonal();
        B.block(0, 0, 2 * N, m - 1) = W.block(0, 0, 2 * N, m - 1);
        B.block(0, m - 1, 2 * N, 2 * N - m) = W.block(0, m, 2 * N, 2 * N - m);
        B.col(2 * N - 1) = -u;
        // compute right hand side
        auto r = W_der * u;
        // solve linear system of equations for derivative of eigenvalue and eigenvector
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu_B(B);
        auto u_der_temp = -lu_B.solve(r);
        // construct eigenvector using derivative of normalization condition
        Eigen::VectorXcd u_der(2 * N);
        u_der.segment(0, m - 1) = u_der_temp.segment(0, m - 1);
        u_der[m - 1] = 0;
        u_der.segment(m, 2 * N - m) = u_der_temp.segment(m - 1, 2 * N - m);
        // compute right hand side for second derivative
        W_der -= (u_der_temp[2 * N - 1] * Eigen::VectorXcd::Ones(2 * N)).asDiagonal();
        auto t = W_der2 * u + 2. * (W_der * u_der);
        // solve linear system of equations second derivative of eigenvector and eigenvalue
        auto u_der2 = -lu_B.solve(t);
        res(2) = u_der2[2 * N - 1].real();
        return res;
    }

}
