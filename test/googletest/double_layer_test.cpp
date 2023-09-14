/**
 * \file double_layer_test.cpp
 * \brief This test file compares Hypersingular BIO for the
 * Helmholtz kernel to a precomputed known solution from
 * a file.
 */
#include <complex>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "galerkin_matrix_builder.hpp"
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
//double c_o = 1;
double c_i = 5;


// define boundary and mesh on boundary
double eps = 0.25;
ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
int numpanels = 50;
ParametrizedMesh mesh(curve.split(numpanels));

// define order of quadrature rule with which to compute matrix entries of operator
unsigned order = 11;

GalerkinMatrixBuilder builder(mesh, discont_space, cont_space, order);
Eigen::VectorXcd K;

// set variables for reading operator from file
Eigen::VectorXcd K_expected(numpanels);
std::ifstream fp_data;
double real, imag;
char sign;
int i = 0;
std::string path = "/home/diego/Uni/Thesis/HelmholtzBEM/raw_data/double_layer_i_" + std::to_string(numpanels) + ".dat";
// use this fie with c_o
//std::string path = "/home/diego/Uni/Thesis/HelmholtzBEM/raw_data/double_layer_o_" + std::to_string(numpanels) + ".dat";

TEST(DoubleLayerTest, compare_row) {
    // compute operator and extract first row
    builder.assembleDoubleLayer(k, c_i);
    K = builder.getDoubleLayer().block(0,0,1,numpanels).transpose();
    // read first row of operator from file
    fp_data.open(path);
    while(fp_data >> real >> imag) {
        K_expected(i) = complex_t((sign=='-')?-real:real,imag);
        i++;
        if (i==numpanels) break;
        fp_data >> sign >> sign;
    }
    std::cout << K.transpose() << std::endl;
    std::cout << "********************" << std::endl;
    std::cout << K_expected.transpose() << std::endl;

    // compare known solution operator from file with computed one
    ASSERT_TRUE((K-K_expected).lpNorm<2>() < 5*1e-4);
    fp_data.close();
}
