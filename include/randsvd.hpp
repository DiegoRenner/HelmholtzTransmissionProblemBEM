/**
 * \file randsvd.hpp
 * \brief This file declares the functions for creating random matrix
 * with standard Gaussian distribution and for approximating the smallest
 * singular value of an invertible matrix by using the randomized SVD
 * approach with subspace iteration for improved accuracy.
 *
 * This file is a part of the HelmholtzTransmissionProblemBEM library.
 *
 * (c) 2023 Luka Marohnić
 */

#ifndef RANDSVD_HPP
#define RANDSVD_HPP

#include <Eigen/Dense>

namespace randomized_svd {

    /**
    * This function creates a random Gaussian matrix with the specified
    * number of rows and columns. The distribution is standard normal N(0,1).
    *
    * @param nr number of rows
    * @param nc number of columns
    */
    Eigen::MatrixXcd randGaussian(int nr, int nc);

    /**
     * This function returns an approximation of the smallest singular value
     * of the given matrix with by using randomized SVD with the given thin
     * random matrix and optionally a number of subspace iterations for better
     * accuracy. LU decomposition with partial pivoting and \f$2q-1\f$
     * Householder QR decompositions are computed, where q is the number of
     * subspace iterations (by default 0).
     *
     * @param T an invertible matrix
     * @param R a thin Gaussian random matrix (2 columns is enough)
     * @param q a positive integer (optional number of subspace iterations)
     */
    double sv(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &R, int q = 0);

    /**
     * This function returns a vector with two elements, an approximation of
     * the smallest singular value of the given matrix and its derivative.
     *
     * @param T an invertible matrix
     * @param T_der the first derivative of T
     * @param R a thin Gaussian random matrix (2 columns is enough)
     * @param q a positive integer (optional number of subspace iterations)
     */
    Eigen::Vector2d sv_der(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &T_der, const Eigen::MatrixXcd &R, int q = 0);

    /**
     * This function returns a vector with three elements, an approximation
     * of the smallest singular value of the given matrix and its first two
     * derivatives.
     *
     * @param T an invertible matrix
     * @param T_der the first derivative of T
     * @param T_der2 the second derivative of T
     * @param R a thin Gaussian random matrix (2 columns is enough)
     * @param q a positive integer (optional number of subspace iterations)
     */
    Eigen::Vector3d sv_der2(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &T_der, const Eigen::MatrixXcd T_der2, const Eigen::MatrixXcd &R, int q);
}

#endif // RANDSVD_HPP
