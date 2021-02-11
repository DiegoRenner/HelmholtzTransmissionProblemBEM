/**
 * \file hypersingular_test.cpp
 * \brief This test file compares the Hypersingular BIO for the
 * Helmholtz kernel to a precomputed known solution from
 * a file.
 */

#include <complex>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "hypersingular.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "abstract_bem_space.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"

typedef std::complex<double> complex_t;

Eigen::MatrixXcd T = Eigen::MatrixXcd::Random(3,3);
Eigen::FullPivLU<Eigen::MatrixXcd> LU(T);
Eigen::MatrixXcd T_herm = T.transpose().conjugate();
Eigen::FullPivLU<Eigen::MatrixXcd> LU_herm(T_herm);

TEST(HypersingularTest, compare_row) {

    std::cout << LU.matrixLU() << std::endl;
    std::cout << "***********************************************" << std::endl;
    std::cout << LU_herm.matrixLU() << std::endl;
    std::cout << "***********************************************" << std::endl;
    std::cout << LU.permutationP().inverse()
                * LU.matrixLU().triangularView<Eigen::UnitLower>().toDenseMatrix()
                * LU.matrixLU().triangularView<Eigen::Upper>().toDenseMatrix()
                * LU.permutationQ().inverse()<< std::endl;
    std::cout << "***********************************************" << std::endl;
    std::cout << T << std::endl;
    std::cout << "***********************************************" << std::endl;
    std::cout << LU.permutationQ()
                 * LU.matrixLU().triangularView<Eigen::Upper>().toDenseMatrix().transpose().conjugate()
                 * LU.matrixLU().triangularView<Eigen::UnitLower>().toDenseMatrix().transpose().conjugate()
                 * LU.permutationP() << std::endl;
}
TEST(Eigenvectors, test) {

}

