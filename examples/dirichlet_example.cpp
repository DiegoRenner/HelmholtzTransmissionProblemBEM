#include <complex>
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_mesh.hpp"
#include "solvers.hpp"
#include "gen_sol.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {

    // define radius of circle
    double eps = 0.25;
    Eigen::Vector2d ipt(0.125,0.0);
    double k = 1.0;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    unsigned order = 11;
    unsigned n_runs = 20;
    double numpanels[n_runs];
    numpanels[0] = 50;
    for (int i=1; i<n_runs; i++){
        numpanels[i] = 2*numpanels[i-1];
    }
    auto fund_sol = [&] (double x1, double x2){
        return sol::fund_sol(k,x1,x2,ipt[0],ipt[1]);
    };
    auto fund_sol_N = [&] (double x1, double x2){
        return sol::fund_sol_N(k,x1,x2,ipt[0],ipt[1]);
    };
    DiscontinuousSpace<0> discont_space;
    ContinuousSpace<1> cont_space;
    // Loop over number of panels
    for (unsigned i = 0; i <= n_runs; i++) {
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        Eigen::VectorXcd Tn_dfk = bvp::direct_first_kind::solve_dirichlet(
                mesh, fund_sol, fund_sol_N, order, k);
    }

    return 0;
}


