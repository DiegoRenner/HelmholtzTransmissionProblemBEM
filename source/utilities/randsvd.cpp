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

namespace randomized_svd {

    Eigen::MatrixXcd randGaussian(int nr, int nc) {
        std::random_device rd {};
        std::mt19937 gen { rd() };
        std::uniform_real_distribution<> d { 0, 1 };
        Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(nr, nc);
        for (int i = 0; i < nr; ++i) for (int j = 0; j < nc; ++j) {
            W(i, j) = std::complex<double>(d(gen), d(gen));
        }
        return W * M_SQRT1_2;
    }

    double sv(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &W, int q) {
        int nr = W.rows(), nc = W.cols();
        Eigen::MatrixXcd Q, thinQ;
        thinQ.setIdentity(nr, nc);
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu_decomp(T);
        Q = lu_decomp.solve(W).householderQr().householderQ() * thinQ;
        for (int i = 0; i < q; ++i) {
            Q = lu_decomp.adjoint().solve(Q).householderQr().householderQ() * thinQ;
            Q = lu_decomp.solve(Q).householderQr().householderQ() * thinQ;
        }
        auto svd = lu_decomp.adjoint().solve(Q).bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        return 1.0 / svd.singularValues()(0);
    }

    Eigen::Vector2d sv_der(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &T_der, const Eigen::MatrixXcd &W, int q) {
        unsigned int N = W.rows();
        Eigen::Vector2d res;
        Eigen::MatrixXcd Wi, Wi_der;
        Wi.setZero(2 * N, 2 * N);
        Wi_der.setZero(2 * N, 2 * N);
        Wi.block(0, N, N, N) = T;
        Wi.block(N, 0, N, N) = T.adjoint();
        Wi_der.block(0, N, N, N) = T_der;
        Wi_der.block(N, 0, N, N) = T_der.adjoint();
        double s = sv(T, W, q);
        Eigen::MatrixXcd A = Wi - (s * Eigen::VectorXcd::Ones(2 * N)).asDiagonal().toDenseMatrix();
        auto qr = A.colPivHouseholderQr();
        Eigen::MatrixXcd Q = qr.householderQ(), R = qr.matrixR();
        Eigen::VectorXcd x = Q.col(2 * N - 1);
        x.normalize();
        res(0) = s;
        res(1) = (x.dot(Wi_der * x)).real();
        return res;
    }

}
