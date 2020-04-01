//
// Created by diegorenner on 2/20/20.
//
#include "hypersingular_smooth_test.hpp"
#include <complex>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "hypersingular_smooth.hpp"
#include "hypersingular_correction.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "abstract_bem_space.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"

typedef std::complex<double> complex_t;
double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

parametricbem2d::DiscontinuousSpace<0> discont_space;
parametricbem2d::ContinuousSpace<1> cont_space;

double k = 1.0;
double n_i = 23.0;
double eps = 0.25;
parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
unsigned order = 11;
double c_o = k;
double c_i = k*sqrt(n_i);
int numpanels = 400;
parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels));
Eigen::VectorXcd W_smooth = parametricbem2d::hypersingular_smooth_helmholtz::GalerkinMatrix(mesh, cont_space, order, c_o).block(0, 0, 1,
                                                                                                                                  numpanels).transpose();
Eigen::VectorXcd W_correction = parametricbem2d::hypersingular_correction::GalerkinMatrix(mesh, cont_space, order).block(0, 0, 1,
                                                                                                                             numpanels).transpose();
Eigen::VectorXcd W = W_smooth + W_correction;
Eigen::VectorXcd W_expected(numpanels);
std::ifstream fp_data;
double real, imag;
char sign;
int i = 0;
std::string path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/raw_data/hypersingular_i_" + std::to_string(numpanels) + ".dat";
TEST(HypersingularTest, disjoint_fair) {
    fp_data.open(path);
    while(fp_data >> real >> imag) {
        W_expected(i) = complex_t((sign=='-')?-real:real,imag);
        i++;
        if (i==numpanels) break;
        fp_data >> sign >> sign;
    }
    fp_data.close();
    std::cout << W.transpose() << std::endl;
    std::cout << "***********************************" << std::endl;
    std::cout << W_expected.transpose() << std::endl;
    std::cout << "***********************************" << std::endl;
    std::cout << (W.segment(2,numpanels-3)-W_expected.segment(2,numpanels-3)).lpNorm<2>()/(numpanels-3) << std::endl;
    ASSERT_TRUE((W.segment(2,numpanels-3)-W_expected.segment(2,numpanels-3)).lpNorm<2>()/(numpanels-3) < sqrt_epsilon);
}

TEST(HypersingularTest, coinciding_precise) {
    fp_data.open(path);
    i = 0;
    while(fp_data >> real >> sign >> imag >> sign) {
        W_expected[i] = complex_t(real, imag);
        i++;
    }
    ASSERT_TRUE(abs(W[0]-W_expected[0]) < 0.0001);
    fp_data.close();
}

TEST(HypersingularTest, adjacent_precise) {
    fp_data.open(path);
    i = 0;
    while(fp_data >> real >> sign >> imag >> sign) {
        W_expected[i] = complex_t(real, imag);
        i++;
    }
    ASSERT_TRUE(abs(W[1]-W_expected[1]) < 0.0001);
    ASSERT_TRUE(abs(W[numpanels-1]-W_expected[numpanels-1]) < 0.0001);
    fp_data.close();
}