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
        std::normal_distribution<> d { 0, 1 };
        Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(nr, nc);
        for (int i = 0; i < nr; ++i) for (int j = 0; j < nc; ++j) {
            W(i, j) = std::complex<double>(d(gen), d(gen));
        }
        return W;
    }

    double sv(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &W, int q) {
        int nr = W.rows(), nc = W.cols();
        Eigen::MatrixXcd Q, B, C, thinQ = Eigen::MatrixXcd::Identity(nr, nc);
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu_decomp(T);
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(lu_decomp.solve(W));
        Q = qr.householderQ() * thinQ;
        for (int i = 0; i < q; ++i) {
            Eigen::HouseholderQR<Eigen::MatrixXcd> qr1(lu_decomp.adjoint().solve(Q));
            Q = qr1.householderQ() * thinQ;
            Eigen::HouseholderQR<Eigen::MatrixXcd> qr2(lu_decomp.solve(Q));
            Q = qr2.householderQ() * thinQ;
        }
        C = Q.transpose().conjugate() * lu_decomp.solve(Q);
        Eigen::BDCSVD<Eigen::MatrixXcd> svd(C * Q.transpose().conjugate());
        return 1.0 / svd.singularValues()(0);
    }

}
