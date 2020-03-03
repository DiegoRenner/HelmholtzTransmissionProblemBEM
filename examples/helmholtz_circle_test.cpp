#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
#include "parametrized_circular_arc.hpp"
#include "parametrized_mesh.hpp"
#include "/home/diegorenner/CLionProjects/Code/third_party/Betl2/Library/analytical_functions/hankel.hpp"
#include "transmission.hpp"
#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"

#define _USE_MATH_DEFINES // for Pi
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {
    double k = 0.9;//679178324354;
    //cout.precision(4);
    //complex_t t = complex_t(1.,1.);
    //std::cout << boost::math::sph_bessel(1,1.2345) << std::endl;
    //std::cout << j1(1.2345) << std::endl;
    double n_o = 1.;
    double n_i = 100.;
    double A_N = 1.;
    double nu = 10;
    //double phi = M_PI/6.;



    complex_t F_nu;
    F_nu = A_N * sqrt(n_i) * (betl2::analytical::sph_bessel_der(nu,k*sqrt(n_i))) *
            (betl2::analytical::sph_hankel1(nu,k))-
           betl2::analytical::sph_hankel1_der(nu,k) *
           boost::math::sph_bessel(nu,k*sqrt(n_i));
    double c_nu = 1/sqrt(M_PI*(boost::math::sph_bessel(nu,k)*boost::math::sph_bessel(nu,k)
            -boost::math::sph_bessel(nu+1,k)*boost::math::sph_bessel(nu-1,k)));
    double c1 = c_nu/(k*k*(n_i-1));
    complex_t c2 = (A_N*(betl2::analytical::sph_bessel_der(nu,k))*
            (betl2::analytical::sph_hankel1(nu,k)) -
            (betl2::analytical::sph_hankel1_der(nu,k))
            *boost::math::sph_bessel(nu,k))/F_nu;

    complex_t c3 = c1/(betl2::analytical::sph_hankel1(nu,k))*
            (boost::math::sph_bessel(nu,k)-boost::math::sph_bessel(nu,k*sqrt(n_i))*c2);



    auto u_inc = [&] (double x1, double x2) {
        double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        return 1./(k*k*(n_o-n_i))*c_nu*boost::math::sph_bessel(nu,k*r)*exp(ii*nu*eta);
    };
    auto u_inc_Tn = [&] (double x1, double x2) {
        double r = sqrt(x1*x1+x2*x2);
        double eta = acos(x1/r);
        Eigen::Vector2cd grad;
        Eigen::Vector2d norm_r;
        Eigen::Vector2d norm_eta;
        Eigen::Vector2d normal;
        // Outward normal vector assuming anticlockwise curve
        normal << r*cos(eta), r*sin(eta);
        // Normalizing the normal vector
        normal /= normal.norm();
        norm_r << cos(eta), sin(eta);
        norm_eta << -sin(eta), cos(eta);
        grad = c_nu*exp(ii*nu*eta)*(k*(betl2::analytical::sph_bessel_der(nu,k*r))*norm_r
                                    +ii*nu/r*boost::math::sph_bessel(nu,k*r)*norm_eta);
        return 1./(k*k*(n_o-n_i))*grad.dot(normal);
    };
    auto u = [&] (double x1, double x2) {
        double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        //complex_t(cos(nu*eta),sin(nu*eta))/(k*k*(n_i-1))*
        //Eigen::Vector2cd ang;
        //ang << exp(ii*x1*cos(phi)), exp(ii*x2*sin(phi));
        return c1*exp(ii*nu*eta)*
               (boost::math::sph_bessel(nu,k*r)-
                boost::math::sph_bessel(nu,k*sqrt(n_i)*r)*c2
               );
    };
    auto u_i = [&] (double x1, double x2) {
        return 2.*u(x1,x2);
    };
    auto u_o = [&] (double x1, double x2) {
       double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        return 2.*c3*exp(ii*nu*eta)*(betl2::analytical::sph_hankel1(nu,k*r));
    };
    auto u_i_Tn = [&] (double x1, double x2) {
        double r = sqrt(x1*x1+x2*x2);
        double eta = acos(x1/r);
        Eigen::Vector2cd grad;
        Eigen::Vector2d norm_r;
        Eigen::Vector2d norm_eta;
        Eigen::Vector2d normal;
        // Outward normal vector assuming anticlockwise curve
        normal << r*cos(eta), r*sin(eta);
        // Normalizing the normal vector
        normal /= normal.norm();
        norm_r << cos(eta), sin(eta);
        norm_eta << -sin(eta), cos(eta);
        grad = 2.*c1*exp(ii*nu*eta)*
                (k*(betl2::analytical::sph_bessel_der(nu,k*r)-
                sqrt(n_i)*(betl2::analytical::sph_bessel_der(nu,k*r*sqrt(n_i)))*c2)*norm_r +
                ii*nu/r*(boost::math::sph_bessel(nu,k*r)-
                        boost::math::sph_bessel(nu,k*r*sqrt(n_i))*c2)*norm_eta);
        return grad.dot(normal);
    };
    auto u_o_Tn = [&] (double x1, double x2) {
        double r = sqrt(x1*x1+x2*x2);
        double eta = acos(x1/r);
        Eigen::Vector2cd grad;
        Eigen::Vector2d norm_r;
        Eigen::Vector2d norm_eta;
        Eigen::Vector2d normal;
        // Outward normal vector assuming anticlockwise curve
        normal << r*cos(eta), r*sin(eta);
        // Normalizing the normal vector
        normal /= normal.norm();
        norm_r << cos(eta), sin(eta);
        norm_eta << -sin(eta), cos(eta);
        grad = -2.*c3* exp(ii*nu*eta)*
                (k*betl2::analytical::sph_hankel1_der(nu,k*r)*norm_r +
                ii*nu/r*betl2::analytical::sph_hankel1(nu,k*r)*norm_eta);
        return grad.dot(normal);
    };
    auto u_inc_T_old = [&] (double x1, double x2) {
        // Getting the tangent vector to evaluate normal
        double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        Eigen::Vector2d normal;
        // Outward normal vector assuming anticlockwise curve
        normal << r*cos(eta), r*sin(eta);
        // Normalizing the normal vector
        normal /= normal.norm();
        Eigen::Vector2cd res;
        complex_t res1;
        complex_t res2;
        double h = 0.1;
        complex_t temp1;
        complex_t temp2;
        Eigen::Vector2d x;
        x << x1, x2;
        temp1 = u_inc(x1-h,x2);
        temp2 = u_inc(x1+h,x2);
        res1 = (temp2-temp1)/(2*h);
        temp1 = u_inc(x1,x2-h);
        temp2 = u_inc(x1,x2+h);
        res2 = (temp2-temp1)/(2*h);
        res << res1, res2;

        return res.dot(normal);
    };
    auto u_inc_T = [&] (double x1, double x2) {
        // Getting the tangent vector to evaluate normal
        double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        Eigen::Vector2d normal;
        // Outward normal vector assuming anticlockwise curve
        normal << r*cos(eta), r*sin(eta);
        // Normalizing the normal vector
        normal /= normal.norm();
        Eigen::Vector2cd res;
        complex_t res1;
        complex_t res2;
        double h = 0.0001;
        double r1;
        double r2;
        complex_t temp1;
        complex_t temp2;
        Eigen::Vector2d x;
        x << x1, x2;
        r1 = sqrt((x1-h)*(x1-h)+x2*x2);
        r2 = sqrt((x1+h)*(x1+h)+x2*x2);
        temp1 = r1 < 1. ? u_inc(x1-h,x2) : complex_t(0.,0.);
        temp2 = r2 < 1. ? u_inc(x1+h,x2) : complex_t(0.,0.);
        res1 = (temp2-temp1)/(2*h);
        r1 = sqrt(x1*x1+(x2-h)*(x2-h));
        r2 = sqrt(x1*x1+(x2+h)*(x2+h));
        temp1 = r1 < 1. ? u_inc(x1,x2-h) : complex_t(0.,0.);
        temp2 = r2 < 1. ? u_inc(x1,x2+h) : complex_t(0.,0.);
        res2 = (temp2-temp1)/(2*h);
        res << res1, res2;
        res1 = (temp2-temp1)/(2*h);

        return res.dot(normal);
    };
    auto u_T = [&] (double x1, double x2) {
        // Getting the tangent vector to evaluate normal
        double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        Eigen::Vector2d normal;
        // Outward normal vector assuming anticlockwise curve
        normal << r*cos(eta), r*sin(eta);
        // Normalizing the normal vector
        normal /= normal.norm();
        Eigen::Vector2cd res;
        complex_t res1;
        complex_t res2;
        double h = 0.1;
        double r1;
        double r2;
        complex_t temp1;
        complex_t temp2;
        Eigen::Vector2d x;
        x << x1, x2;
        r1 = sqrt((x1-h)*(x1-h)+x2*x2);
        r2 = sqrt((x1+h)*(x1+h)+x2*x2);
        temp1 = r1 < 1. ? u_i(x1-h,x2) : u_o(x1-h,x2);
        temp2 = r2 < 1. ? u_i(x1+h,x2) : u_o(x1+h,x2);
        res1 = (temp2-temp1)/(2*h);
        r1 = sqrt(x1*x1+(x2-h)*(x2-h));
        r2 = sqrt(x1*x1+(x2+h)*(x2+h));
        temp1 = r1 < 1. ? u_i(x1,x2-h) : u_o(x1,x2-h);
        temp2 = r2 < 1. ? u_i(x1,x2+h) : u_o(x1,x2+h);
        res2 = (temp2-temp1)/(2*h);
        res << res1, res2;
        res1 = (temp2-temp1)/(2*h);

        return res.dot(normal);
    };
    auto sol_d = [&] (double x1, double x2) {
        return u_i(x1,x2)-u_inc(x1,x2);
    };
    auto sol_n = [&] (double x1, double x2) {
        return u_i_Tn(x1,x2)-u_inc_Tn(x1,x2);
    };
/*
    double alpha = 5/3*M_PI;
    double x1 = cos(alpha);
    double x2 = sin(alpha);
    // std::cout << u(x1,x2) << " " << c1 << " " << c2 << " " << F_nu << std::endl;
    // Plot u_inc
    int n = 201;
    Eigen::MatrixXd u_plot(n,n);
    double x;
    double y;
    double r;
    double eta;
    for (int j = 0; j < n; ++j) {
        for (int l = 0; l < n; ++l) {
            x = 2*M_PI*double(j-n/2)/(n-1);
            y = 2*M_PI* double(l-n/2)/(n-1);
            r = sqrt(x*x+y*y);
            eta = (r==0.0) ? 0 : acos(x/r);
            /*std::cout <<"x: " << x << std::endl;
            std::cout <<"y: " << y << std::endl;
            std::cout <<"r: " << r << std::endl;
            std::cout <<"eta: " << eta << std::endl;
            if ( r < 1.) {
                u_plot(j, l) = std::abs(u_inc(x,y));
            } else {
                u_plot(j, l) = 0.;
            }
            //u_inc(j,l) = std::abs(u_inc(r,eta));
        }
    }
    //std::cout << u_plot;*/

    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),0.5 ,0,2*M_PI);
    unsigned maxpanels = 200;
    unsigned order = 11;
    double c_o = k*sqrt(n_o);
    double c_i = k*sqrt(n_i);

    parametricbem2d::DiscontinuousSpace<0> discont_space;
    parametricbem2d::ContinuousSpace<1> cont_space;
    // Loop over number of panels
    double nps[4];
    nps[0] = 50.;
    nps[1] = 100.;
    nps[2] = 200.;
    nps[3] = 500.;
    for (unsigned numpanels = 0; numpanels <= 3; numpanels += 1) {
        parametricbem2d::ParametrizedMesh mesh(curve.split(nps[numpanels]));
        Eigen::VectorXcd Tn_dfk = parametricbem2d::transmission_bvp::direct_second_kind::solve(
                mesh, u_inc, u_inc_Tn, u_i, u_i_Tn, order, c_o, c_i);
        /*Eigen::VectorXcd u_exp_dir = cont_space.Interpolate_helmholtz(u_o, mesh);
        Eigen::VectorXcd u_exp_n = discont_space.Interpolate_helmholtz(u_o_Tn, mesh);
        Eigen::VectorXcd  u_exp_N(2*numpanels);
        u_exp_N << u_exp_dir, u_exp_n;
        Eigen::VectorXcd res(2*numpanels);
        using PanelVector = parametricbem2d::PanelVector;
        PanelVector panels = mesh.getPanels();
        for (int i = 0; i < numpanels; i++) {
            Eigen::Vector2d pt = panels[i]->operator()(0.);
            res(i) = (Tn_dfk(i) - u_o(pt(0),pt(1)));
            res(i+numpanels) = (Tn_dfk(i+numpanels) - u_o_Tn(pt(0),pt(1)) );
        }*/
        //std::ofstream filename;
        //filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_data.dat", std::ios_base::app);
        //filename <<  res.segment(0,numpanels).lpNorm<Eigen::Infinity>() << " " << sqrt(2*M_PI/numpanels)
        //         << " " << res.segment(0,numpanels).lpNorm<Eigen::Infinity>()/sqrt(2*M_PI/numpanels) << " " <<
        //         u_exp_dir.cwiseAbs().mean() << std::endl;

        //std::ofstream filename1;
        //filename1.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_data.dat", std::ios_base::app);
        //filename1 <<  res.segment(numpanels,numpanels).lpNorm<Eigen::Infinity>() << " " << sqrt(2*M_PI/numpanels)
        //         << " " << res.segment(numpanels,numpanels).lpNorm<Eigen::Infinity>()/sqrt(2*M_PI/numpanels) << " " <<
        //         u_exp_n.cwiseAbs().mean() << std::endl;

        //std::cout << res.lpNorm<Eigen::Infinity>() << std::endl;
        //std::cout << u_exp_dir<< std::endl;
        //std::cout << "-------" << std::endl;
        //std::cout << u_exp_n<< std::endl;
        //std::cout << "-------" << std::endl;
        //std::cout << Tn_dfk.segment(0,numpanels) << std::endl;
        //std::cout << "-------" << std::endl;
        //std::cout << Tn_dfk.segment(numpanels,numpanels) << std::endl;
        //std::cout << "-------" << std::endl;
    }


    // Verifying that Galerkin matrix is circulant
/*  for (unsigned int i = 0; i < N - 1; ++i) {
    Eigen::VectorXcd col1 = galerkin.col(i);
    Eigen::VectorXcd col2(N);
    col2 << galerkin.col(i + 1).segment(1, N - 1), galerkin.col(i + 1)(0);
    assert((col1 - col2).norm() < 1e-7);
  }
  // Getting the eigenvalues
  Eigen::EigenSolver<Eigen::MatrixXcd> es(galerkin);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> eigs =
      es.eigenvalues();
  Eigen::VectorXd eigs_r(N);
  // Extracting the real part of eigenvalues (imaginary zero in this case)
  for (unsigned int i = 0; i < N; ++i)
    eigs_r(i) = eigs(i).real();
  // Sorting the real parts of eigenvalues
  std::sort(eigs_r.data(), eigs_r.data() + eigs_r.size());
  // Saving the real parts to a file
  for (unsigned int i = 0; i < N; ++i)
    output << eigs_r(i) << std::endl;*/
    return 0;
}


