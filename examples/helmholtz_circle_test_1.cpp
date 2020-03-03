//#include <cmath>
#include <complex>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <math.h>
//
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
#include "parametrized_circular_arc.hpp"
#include "transmission.hpp"
#include "generate_solution.hpp"

#define _USE_MATH_DEFINES // for Pi
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {

    // define radius of circle
    double eps = 0.25;
    int l = 0;
    double a_n[2*l+1];
    Eigen::Vector2d ipt(0.125,0.0);
    //a_n[1] = 1.0;
    //a_n[2] = 1.0;
    double k = 0.63;
    double n_i = 23.;
    double n_o = 1.;
    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    unsigned order = 11;
    double c_o = k*sqrt(n_o);
    double c_i = k*sqrt(n_i);
    a_n[0] = 1./((k*k*(n_o-n_i))*sqrt(M_PI*eps*eps*(jn(0,k)*jn(0,k)-jn(-1,k)*jn(1,k))));
    unsigned n_runs = 10;
    double numpanels[n_runs];
    for (int i=0; i<n_runs; i++){
        numpanels[i] = (i+11)*50;
    }

    auto f_i = [&] (double x1, double x2) {
        return sol::u_i(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_s = [&] (double x1, double x2) {
        return sol::u_s(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_t = [&] (double x1, double x2) {
        return sol::u_t(x1, x2, l, eps, a_n, k, n_i);
    }
;

    auto f_i_N = [&] (double x1, double x2) {
        return sol::u_i_N(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_s_N = [&] (double x1, double x2) {
        return sol::u_s_N(x1, x2, l, eps, a_n, k, n_i);
    };

    auto u_t_N = [&] (double x1, double x2) {
        return sol::u_t_N(x1, x2, l, eps, a_n, k, n_i);
    };

    auto sol = [&] (double x1, double x2){
            return (sol::u_i(x1, x2,l,eps,a_n,k,n_i)
                    + sol::u_s(x1, x2, l, eps, a_n, k, n_i));
    };
    auto sol_N = [&] (double x1, double x2){
            return (sol::u_i_N(x1,x2,l,eps,a_n,k,n_i)
                  + sol::u_s_N(x1, x2, l, eps, a_n, k, n_i));
    };

    auto fund_sol = [&] (double x1, double x2){
        return sol::fund_sol(k,x1,x2,ipt[0],ipt[1]);
    };

    auto fund_sol_N = [&] (double x1, double x2){
        return sol::fund_sol_N(k,x1,x2,ipt[0],ipt[1]);
    };
    parametricbem2d::DiscontinuousSpace<0> discont_space;
    parametricbem2d::ContinuousSpace<1> cont_space;
    // Loop over number of panels
    for (unsigned i = 0; i <= n_runs; i++) {
        parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels[i]));
        Eigen::VectorXcd Tn_dfk = parametricbem2d::transmission_bvp::direct_second_kind::solve(
                mesh, f_i, f_i_N, u_t, u_t_N, order, c_o, c_i);
    }

    return 0;
}


