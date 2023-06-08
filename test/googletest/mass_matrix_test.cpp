/**
 * \file mass_matrix_test.cpp
 * \brief This test file compares the Identitiy BIO
 * for different FEM-spaces to a precomputed known
 * solution from a file.
 */
#include "mass_matrix.hpp"
#include <complex>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "abstract_bem_space.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"

typedef std::complex<double> complex_t;

// define FEM-spaces of lowest order
DiscontinuousSpace<0> discont_space;
ContinuousSpace<1> cont_space;

// define wavenumber and refraction index
double k = 1.0;
double n_i = 23.0;
double c_o = k;
double c_i = k*sqrt(n_i);

// define boundary and mesh on boundary
double eps = 0.25;
ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
int numpanels = 50;
ParametrizedMesh mesh(curve.split(numpanels));

// define order of quadrature rule with which to compute matrix entries of operator
unsigned order = 11;

// run test for continuous FEM-sapces
TEST(SingleLayerTest, cont_sapce) {
    // compute operator and extract first row
    Eigen::VectorXcd M = mass_matrix::GalerkinMatrix(mesh, cont_space, cont_space, getGaussQR(order, 0., 1.)).block(0,0,1,numpanels).transpose();

    // set variables for reading operator from file
    Eigen::VectorXcd M_expected(numpanels);
    std::ifstream fp_data;
    double real;
    char sign;
    int i = 0;
    std::string path =
            "/home/diegorenner/Uni/Thesis/HelmholtzBEM/raw_data/mass_matrix_cont_" + std::to_string(numpanels) + ".dat";
    // read first row of operator from file
    fp_data.open(path);
    while(fp_data >> real) {
        M_expected(i) = real;
        i++;
        if (i==numpanels) break;
        fp_data >> sign;
    }

    // compare known solution operator from file with computed one
    ASSERT_TRUE((M-M_expected).lpNorm<2>() < 5e-5);
    fp_data.close();
}

// run test for discontinuous FEM-sapces
TEST(SingleLayerTest, discont_spaces) {
    // compute operator and extract first row
    Eigen::VectorXcd M = mass_matrix::GalerkinMatrix(mesh, cont_space, discont_space, getGaussQR(order, 0., 1.)).block(0,0,1,numpanels).transpose();

    // set variables for reading operator from file
    Eigen::VectorXcd M_expected(numpanels);
    std::ifstream fp_data;
    double real;
    char sign;
    int i = 0;
    std::string path =
            "/home/diegorenner/Uni/Thesis/HelmholtzBEM/raw_data/mass_matrix_discont_" + std::to_string(numpanels) + ".dat";
    // read first row of operator from file
    fp_data.open(path);
    while(fp_data >> real) {
        M_expected(i) = real;
        i++;
        if (i==numpanels) break;
        fp_data >> sign;
    }


    // compare known solution operator from file with computed one
    std::cout << (M-M_expected).lpNorm<2>() << std::endl;
    ASSERT_TRUE((M-M_expected).lpNorm<2>() < 5e-5);
    fp_data.close();
}
