//
// Created by diegorenner on 2/20/20.
//
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
double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

DiscontinuousSpace<0> discont_space;
ContinuousSpace<1> cont_space;

double k = 1.0;
double n_i = 23.0;
double eps = 0.25;
ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
unsigned order = 11;
double c_o = k;
double c_i = k*sqrt(n_i);
int numpanels = 50;
ParametrizedMesh mesh(curve.split(numpanels));
Eigen::VectorXcd V = single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k, n_i).block(0,0,1,numpanels).transpose();
Eigen::VectorXcd V_expected(numpanels);
std::ifstream fp_data;
double real, imag;
char sign;
int i = 0;
std::string path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/raw_data/single_layer_i_" + std::to_string(numpanels) + ".dat";

TEST(SingleLayerTest, disjoint_fair) {

    fp_data.open(path);
    while(fp_data >> real >> imag) {
        V_expected(i) = complex_t((sign=='-')?-real:real,imag);
        i++;
        if (i==numpanels) break;
        fp_data >> sign >> sign;
    }
    std::cout << V.transpose() << std::endl;
    std::cout << "***********************************" << std::endl;
    std::cout << V_expected.transpose() << std::endl;
    std::cout << "***********************************" << std::endl;
    std::cout << (V.segment(2,numpanels-3)-V_expected.segment(2,numpanels-3)).lpNorm<2>()/(numpanels-3) << std::endl;
    ASSERT_TRUE((V.segment(2,numpanels-3)-V_expected.segment(2,numpanels-3)).lpNorm<2>()/(numpanels-3) < sqrt_epsilon);
    //fp_data.close();
}

TEST(SingleLayerTest, coinciding_precise) {
    fp_data.open(path);
    i = 0;
    while(fp_data >> real >> sign >> imag >> sign) {
        V_expected[i] = complex_t(real, imag);
        i++;
    }
    ASSERT_TRUE(abs(V[0]-V_expected[0]) < 0.0001);
    fp_data.close();
}

TEST(SingleLayerTest, adjacent_precise) {
    fp_data.open(path);
    i = 0;
    while(fp_data >> real >> sign >> imag >> sign) {
        V_expected[i] = complex_t(real, imag);
        i++;
    }
    ASSERT_TRUE(abs(V[1]-V_expected[1]) < 0.0001);
    ASSERT_TRUE(abs(V[numpanels-1]-V_expected[numpanels-1]) < 0.0001);
    fp_data.close();
}
