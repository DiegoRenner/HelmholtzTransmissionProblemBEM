#include <complex>
#include <iostream>
#include "continuous_space.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_mesh.hpp"
#include "solvers.hpp"
#include "gen_sol.hpp"
#include "mass_matrix.hpp"

/**
 * This script computes the solution to a Neumann problem
 * using first kind direct BIEs.
 * No command line parameters are necessary.
 */
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {

    // define boundary
    double eps = 0.25;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    //defin wavenumber
    double k = 0.63;
    // define order of quadrature rule with which to to compute BIOs
    unsigned order = 11;
    // define #panels with which to compute BIOs
    unsigned n_runs = 7;
    double numpanels[n_runs];
    numpanels[0] = 50;
    for (int i=1; i<n_runs; i++){
        numpanels[i] = 2*numpanels[i-1];
    }
    // define functions for Dirichlet and Neumann data
    Eigen::Vector2d ipt(0.125,0.0);
    auto fund_sol_dir = [&] (double x1, double x2){
        return sol::fund_sol_dir(k, x1, x2, ipt[0], ipt[1]);
    };
    auto fund_sol_neu = [&] (double x1, double x2){
        return sol::fund_sol_neu(k, x1, x2, ipt[0], ipt[1]);
    };
    // define FEM-sapces for result validation later on
    ContinuousSpace<1> cont_space;
    // Loop over #panels
    for (unsigned i = 0; i <= n_runs; i++) {
        // compute mesh
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        // compute Dirichlet data from Neumann data using 1st kind BIE
        Eigen::VectorXcd res = bvp::direct_first_kind::solve_neumann(
                mesh, fund_sol_neu, order, k);
        // compute interpolation coefficients of known solution Dirichlet data in FEM-sapce
        Eigen::VectorXcd res_known = cont_space.Interpolate_helmholtz(fund_sol_dir,mesh);
        // compute mass matrix for projection onto orthonormal basis functions
        Eigen::MatrixXcd M = mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
        // update user on residual error between computed and known FEM-space interpolation coefficients
        // that have been projected onto orthornomal basis
        std::cout << "Computed Dirichlet data on " << numpanels[i] << " panels." << std::endl;
        std::cout << "Residual error of FEM-space interpolation coefficients:" << std::endl;
        std::cout << sqrt(((res-res_known)).dot(M*(res-res_known))).real() << std::endl;
        std::cout << "*******************************************************" << std::endl;
    }

    return 0;
}


