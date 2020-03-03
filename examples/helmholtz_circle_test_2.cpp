//#include <algorithm>
//#include <cassert>
#include <cmath>
#include <complex>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <math.h>

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
#include "parametrized_circular_arc.hpp"
#include "parametrized_mesh.hpp"
#include "/home/diegorenner/CLionProjects/Code/third_party/Betl2/Library/analytical_functions/hankel.hpp"
#include "transmission.hpp"
//#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "generate_solution.hpp"

#define _USE_MATH_DEFINES // for Pi
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {
    double k = 0.9;//679178324354;
    double n_o = 1.;
    double n_i = 100.;
    double A_N = 1.;
    double nu = 10;



    complex_t F_nu;
    F_nu = A_N * sqrt(n_i) * (sol::jn_der(nu,k*sqrt(n_i))) *(sol::hn(nu,k)) - sol::hn_der(nu,k) * jn(nu,k*sqrt(n_i));
    double c_nu = 1/sqrt(M_PI*(jn(nu,k)*jn(nu,k)-jn(nu+1,k)*jn(nu-1,k)));
    double c1 = c_nu/(k*k*(n_i-1));
    complex_t c2 = (A_N*(sol::jn_der(nu,k))*
                    (sol::hn(nu,k)) - (sol::hn_der(nu,k))*jn(nu,k))/F_nu;
    //complex_t c3 = c1/(betl2::analytical::sph_hankel1(nu,k))*
    //        (boost::math::sph_bessel(nu,k)-boost::math::sph_bessel(nu,k*sqrt(n_i))*c2);

    auto u_inc = [&] (double x1, double x2) {
        double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        return 1./(k*k*(n_o-n_i))*c_nu*jn(nu,k*r)*exp(ii*nu*eta);
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
        grad = c_nu*exp(ii*nu*eta)*(k*(sol::jn_der(nu,k*r))*norm_r
                                    +ii*nu/r*jn(nu,k*r)*norm_eta);
        return 1./(k*k*(n_o-n_i))*grad.dot(normal);
    };
    auto u = [&] (double x1, double x2) {
        double r = sqrt(x1*x1 + x2*x2);
        double eta = acos(x1/r);
        //complex_t(cos(nu*eta),sin(nu*eta))/(k*k*(n_i-1))*
        //Eigen::Vector2cd ang;
        //ang << exp(ii*x1*cos(phi)), exp(ii*x2*sin(phi));
        return c1*exp(ii*nu*eta)*
               (jn(nu,k*r)-
                jn(nu,k*sqrt(n_i)*r)*c2
               );
    };
    auto u_i = [&] (double x1, double x2) {
        return u(x1,x2);
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
        grad = c1*exp(ii*nu*eta)*
                (k*(sol::jn_der(nu,k*r)-
                sqrt(n_i)*(sol::jn_der(nu,k*r*sqrt(n_i)))*c2)*norm_r +
                ii*nu/r*(jn(nu,k*r)-jn(nu,k*r*sqrt(n_i))*c2)*norm_eta);
        return grad.dot(normal);
    };

    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),1.0 ,0,2*M_PI);
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
        Eigen::VectorXcd Tn_dfk = parametricbem2d::transmission_bvp::direct_second_kind::solve_1(
                mesh, u_inc, u_inc_Tn, u_i, u_i_Tn, order, c_o, c_i);
    }

    return 0;
}


