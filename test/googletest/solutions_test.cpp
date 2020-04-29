//
// Created by diegorenner on 2/20/20.
//
#include <complex>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "double_layer.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>
//#include <CLucene/search/Query.h>
#include "parametrized_circular_arc.hpp"
#include "abstract_bem_space.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "generate_solution.hpp"
#include "math.h"
#include <boost/math/special_functions.hpp>
#include </usr/include/complex_bessel.h>

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0.,1.);
double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
double tol = std::numeric_limits<double>::epsilon();

parametricbem2d::DiscontinuousSpace<0> discont_space;
parametricbem2d::ContinuousSpace<1> cont_space;

double k = 0.63;
double n_i = 23.0;
double eps = 0.7;
int l = 0;
double a_n[1];
parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
double c_o = k;
double c_i = k*sqrt(n_i);
int numpanels = 10;
parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels));

TEST(FuntionsTest, Bessel){
    ASSERT_TRUE(jn(0,0.0)==1.0);
    ASSERT_TRUE(boost::math::cyl_bessel_j(0,0.0)==1.0);
    ASSERT_TRUE(sp_bessel::besselJ(0,0.0)==1.0);

    ASSERT_NEAR(jn(0,M_PI),-0.304242177644093864202,tol);
    ASSERT_NEAR(boost::math::cyl_bessel_j(0,M_PI),-0.304242177644093864202,tol);
    ASSERT_NEAR(sp_bessel::besselJ(0,M_PI,0.0).real(),-0.304242177644093864202,tol);

    ASSERT_NEAR(sol::jn_der(0.0,0.0).real(),0.0,tol);
    ASSERT_NEAR(boost::math::cyl_bessel_j_prime(0.0,0.0),0.0,tol);
    ASSERT_NEAR(sp_bessel::besselJp(0.0,0.0).real(),0.0,tol);

    ASSERT_NEAR(sol::jn_der(0.0,M_PI).real(),-0.284615343179752757345310599686,tol);
    ASSERT_NEAR(boost::math::cyl_bessel_j_prime(0.0,M_PI),-0.284615343179752757345310599686,tol);
    ASSERT_NEAR(sp_bessel::besselJp(0.0,M_PI).real(),-0.284615343179752757345310599686,tol);

    ASSERT_NEAR(sol::hn_der(0.0,0.1).imag(),6.45895,100*sqrt_epsilon);
    ASSERT_NEAR(sol::hn_der(0.0,0.1).real(),-0.0499375,100*sqrt_epsilon);
    ASSERT_NEAR(boost::math::cyl_neumann_prime(0.0,0.1),6.45895,100*sqrt_epsilon);
    ASSERT_NEAR(boost::math::cyl_bessel_j_prime(0.0,0.1),-0.0499375,100*sqrt_epsilon);
    ASSERT_NEAR(sp_bessel::hankelH1p(0.0,0.1).imag(),6.45895,100*sqrt_epsilon);
    ASSERT_NEAR(sp_bessel::hankelH1p(0.0,0.1).real(),-0.0499375,100*sqrt_epsilon);

    ASSERT_NEAR(exp(ii*(3./4.)*M_PI).imag(),complex_t(cos((3./4.)*M_PI),sin((3./4.)*M_PI)).imag(),tol);
    ASSERT_NEAR(sol::r_coeff(30,1.0,31.4683226,4.).real(),-1.,sqrt_epsilon);

    //ASSERT_NEAR(sol::jn_der(0.0,M_PI),-0.284615343179752757345310599686,tol);
};

TEST(SolutionsTest, dirichlet) {
    //for ( int i = 0; i < 2*l+1; i++) {
    //    a_n[i] = 0.0;
    //}
    //a_n[60] = 1.0;
    a_n[0] = 1.0;
    //a_n[1] = 1.0;
    //a_n[2] = 1.0;

    auto u_i = [&] (double x1, double x2) {
        return sol::u_i(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_s = [&] (double x1, double x2) {
        return sol::u_s(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_t = [&] (double x1, double x2) {
        return sol::u_t(x1, x2, l, eps, a_n, k, n_i);
    };

    Eigen::VectorXcd outer = cont_space.Interpolate_helmholtz(u_s, mesh);
    Eigen::VectorXcd inner = cont_space.Interpolate_helmholtz(u_t, mesh);
    Eigen::VectorXcd incoming = cont_space.Interpolate_helmholtz(u_i, mesh);
    std::cout << inner << std::endl;
    std::cout << "***************************************" << std::endl;
    std::cout << (outer + incoming) << std::endl;
    std::cout << "***************************************" << std::endl;
    std::cout << u_i(0.5,0.0) << std::endl;
    std::cout << u_s(0.5,0.0) + u_i(0.5,0.0)<< std::endl;
    std::cout << u_t(0.5,0.0) << std::endl;
    ASSERT_TRUE(inner.isApprox((outer+incoming)));
}

TEST(SolutionsTest, neumann) {
    //for ( int i = 0; i < 2*l+1; i++) {
    //    a_n[i] = 0.0;
    //}
    //a_n[60] = 1.0;
    a_n[0] = 1.0;
    //a_n[1] = 1.0;
    //a_n[2] = 1.0;

    auto u_i_N = [&] (double x1, double x2) {
        return sol::u_i_N(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_s_N = [&] (double x1, double x2) {
        return sol::u_s_N(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_t_N = [&] (double x1, double x2) {
        return sol::u_t_N(x1, x2, l, eps, a_n, k, n_i);
    };
    Eigen::VectorXcd outer = discont_space.Interpolate_helmholtz(u_s_N, mesh);
    Eigen::VectorXcd inner = discont_space.Interpolate_helmholtz(u_t_N, mesh);
    Eigen::VectorXcd incoming = discont_space.Interpolate_helmholtz(u_i_N, mesh);
    std::cout << inner << std::endl;
    std::cout << "***************************************" << std::endl;
    std::cout << k*eps*(outer + incoming) << std::endl;
    std::cout << "***************************************" << std::endl;
    std::cout << u_i_N(0.5,0.0) << std::endl;
    std::cout << k*eps*(u_s_N(0.5,0.0) + u_i_N(0.5,0.0) )<< std::endl;
    std::cout << u_t_N(0.5,0.0) << std::endl;
    ASSERT_TRUE(inner.isApprox((outer+incoming)));
}
