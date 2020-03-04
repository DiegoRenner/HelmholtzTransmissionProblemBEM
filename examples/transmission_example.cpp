#include <complex>
#include "parametrized_circular_arc.hpp"
#include "solvers.hpp"
#include "generate_solution.hpp"
#include "discontinuous_space.hpp"
#include "continuous_space.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {

    // define radius of circle
    double eps = 0.25;
    int l = 0;
    double a_n[2*l+1];
    Eigen::Vector2d ipt(0.125,0.0);
    double n_i = 23.;
    double n_o = 1.;
    double k = 0.63;
    double k_i = k*sqrt(n_i);
    double k_o = k*sqrt(n_o);
    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    unsigned order = 11;
    a_n[0] = 1./((k*k*(n_o-n_i))*sqrt(M_PI*eps*eps*(jn(0,k)*jn(0,k)-jn(-1,k)*jn(1,k))));
    unsigned n_runs = 20;
    double numpanels[n_runs];
    numpanels[0] = 50;
    for (int i=1; i<n_runs; i++){
        numpanels[i] = 2*numpanels[i-1];
    }
    auto u_i_dir = [&] (double x1, double x2) {
        return sol::u_i(x1, x2, l, eps, a_n, k, n_i);
    };
    auto u_t_dir = [&] (double x1, double x2) {
        return sol::u_t(x1, x2, l, eps, a_n, k, n_i);
    };
    auto u_i_neu = [&] (double x1, double x2) {
        return sol::u_i_N(x1, x2, l, eps, a_n, k, n_i);
    };
    auto u_t_neu = [&] (double x1, double x2) {
        return sol::u_t_N(x1, x2, l, eps, a_n, k, n_i);
    };
    parametricbem2d::DiscontinuousSpace<0> discont_space;
    parametricbem2d::ContinuousSpace<1> cont_space;
    // Loop over number of panels
    for (unsigned i = 0; i <= n_runs; i++) {
        parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels[i]));
        Eigen::VectorXcd Tn_dfk = parametricbem2d::tsp::direct_second_kind::solve(
                mesh, u_i_dir, u_i_neu, u_t_dir, u_t_neu, order, k_o, k_i);
    }
return 0;
}
