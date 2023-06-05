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
        Eigen::MatrixXcd Q, thinQ = Eigen::MatrixXcd::Identity(nr, nc);
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu_decomp(T);
        Q = lu_decomp.solve(W).householderQr().householderQ() * thinQ;
        for (int i = 0; i < q; ++i) {
            Q = lu_decomp.adjoint().solve(Q).householderQr().householderQ() * thinQ;
            Q = lu_decomp.solve(Q).householderQr().householderQ() * thinQ;
        }
        return 1.0 / lu_decomp.adjoint().solve(Q).bdcSvd().singularValues()(0);
    }

}
