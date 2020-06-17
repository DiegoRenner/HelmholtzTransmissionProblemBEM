/**
 * \file dirichlet_example.cpp
 * \brief This target builds a script that computes the solution to a Dirichlet problem
 * using first kind direct BIEs.
 * No command line parameters are necessary.
 * Once built the script can be run as follows:
 *
 * <tt>
 * /path/to/library/bin/dirichlet_example.
 * </tt>
 *
 * The user will be updated over the residual error in 
 * the euclidean norm of the computed FEM-space 
 * interpolation coefficients to the known FEM-space 
 * interpolation coefficients for the current number of panels through the command line.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */
#include <complex>
#include <iostream>
#include "discontinuous_space.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_mesh.hpp"
#include "solvers.hpp"
#include "gen_sol.hpp"
#include "mass_matrix.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {
    // define boundary
    double eps = 0.25;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    // define wavenumber
    double k = 1.0;
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
        return sol::fund_sol_dir(k,x1,x2,ipt[0],ipt[1]);
    };
    auto fund_sol_neu = [&] (double x1, double x2){
        return sol::fund_sol_neu(k,x1,x2,ipt[0],ipt[1]);
    };
    // define FEM-sapces for result validation later on
    DiscontinuousSpace<0> discont_space;
    // Inform user of started computation.
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Solving Dirichlet problem for increasing grid sizes." << std::endl;
    std::cout << "Using first-kind direct BIEs." << std::endl;
    std::cout << "Using lowest order FEM-spaces." << std::endl;
    std::cout << std::endl;
    // loop over #panels
    for (unsigned i = 0; i <= n_runs; i++) {
        // compute mesh
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        // compute Neumann data from Dirichlet data using 1st kind BIE
        Eigen::VectorXcd res = bvp::direct_first_kind::solve_dirichlet(
                mesh, fund_sol_dir, order, k);
        // compute interpolation coefficients of known solution Neumann data in FEM-sapce
        Eigen::VectorXcd res_known = discont_space.Interpolate_helmholtz(fund_sol_neu,mesh);
        // compute mass matrix for projection onto orthonormal basis functions
        Eigen::MatrixXcd M = mass_matrix::GalerkinMatrix(mesh,discont_space,discont_space,order);
        // update user on residual error between computed and known FEM-space 
	// interpolation coefficients that have been projected onto orthornomal basis
        std::cout << "#######################################################" << std::endl;
        std::cout << "Computed Neumann data on " << numpanels[i] << " panels." << std::endl;
        std::cout << "Residual error of FEM-space interpolation coefficients:" << std::endl;
        std::cout << sqrt(((res-res_known)).dot(M*(res-res_known))).real() << std::endl;
        std::cout << "#######################################################" << std::endl;
        std::cout << std::endl;

    }

    return 0;
}


