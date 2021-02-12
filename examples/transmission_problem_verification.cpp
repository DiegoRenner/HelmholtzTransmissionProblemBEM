/**
 * \file transmission_problem_verification.cpp
 * \brief This target builds a script that computes solutions to
 * the analytically solvable case of the Helmholtz transmission
 * problem where the scatterer is a circle using second-kind direct
 * BIEs and Galerkin BEM.
 * The results are written to file.
 * The script can be run as follows:
 * <tt>
 * /path/to/transmission_problem_verification \<radius of circle\>
 *      \<\#coeffs for series expansion of solution\> \<refraction inside\>
 *      \<refraction outside\> \<wavenumber\>
 *      \<order of quadrature rule\> \<outputfile\>
 * </tt>
 * This output file will contain two columns.
 * The first will contain the current panel size.
 * The second will contain the residual error in the euclidean 
 * norm of the computed FEM-space interpolation coefficients to the 
 * known FEM-space interpolation coefficients for the current number of panels.
 * Then the columns will contain the computed derivative, the 
 * extrapolated derivative, the computed second derivative and the extrapolated second 
 * derivative in this order.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "solvers.hpp"
#include "gen_sol.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "mass_matrix.hpp"
#include <chrono>
using namespace std::chrono;

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main(int argc, char** argv) {

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[3]);
    double c_o = atof(argv[4]);
    double k = atof(argv[5]);

    //set boundary
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);

    // define order of quadrature rule used to compute matrix entries
    unsigned order = atoi(argv[6]);

    // clear existing file
    std::ofstream filename;
    filename.open(argv[7], std::ofstream::out | std::ofstream::trunc);
    filename.close();

    // set number of panels for which to compute BIOs
    unsigned n_runs = 7;
    int numpanels[n_runs];
    numpanels[0] = 50;
    for (int i=1; i<n_runs; i++){
        numpanels[i] = 2*numpanels[i-1];
    }
    // define incoming wave and resulting wave
    int l = atoi(argv[2]);
    double a_n[2*l+1];
    for( int i = 0; i<2*l+1; i++) {
        a_n[i] = 1./((k*k*(c_o-c_i))*sqrt((2*l+1)*M_PI*eps*eps*(jn(i-l,k)*jn(i-l,k)-jn(i-l-1,k)*jn(i-l+1,k))));
    }
    auto u_i_dir = [&] (double x1, double x2) {
        return sol::u_i(x1, x2, l, a_n, k);
    };
    auto u_t_dir = [&] (double x1, double x2) {
        return sol::u_t(x1, x2, l, eps, a_n, k, c_i);
    };
    auto u_i_neu = [&] (double x1, double x2) {
        return sol::u_i_neu(x1, x2, l, a_n, k);
    };
    auto u_t_neu = [&] (double x1, double x2) {
        return sol::u_t_neu(x1, x2, l, eps, a_n, k, c_i);
    };

    // set FEM-sapces of lowest order for validation
    ContinuousSpace<1> cont_space;

	#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Computing convergence rate for Helmholtz TP BIEs." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
	#endif

    // Loop over number of panels
    for (unsigned i = 0; i < n_runs; i++) {

        // generate mesh on boudary
        ParametrizedMesh mesh(curve.split(numpanels[i]));

        // compute interpolation coefficients
        // in FEM-sapces for resulting waves
        // DEBUG auto start = high_resolution_clock::now();
        Eigen::VectorXcd sol = tp::direct_second_kind::solve(
                mesh, u_i_dir, u_i_neu, order, k, c_o, c_i);
        // DEBUG auto end = high_resolution_clock::now();
        // DEBUG auto duration = duration_cast<milliseconds>(end - start);

        // compute mass matrix for projection on to orthonromal basis functions
        Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
        Eigen::MatrixXcd M(2*numpanels[i],2*numpanels[i]);
        M.block(0,0,numpanels[i],numpanels[i]) = M_cont;
        M.block(numpanels[i],numpanels[i],numpanels[i],numpanels[i]) = M_cont;

        // compute interpolation coefficients in FEM-spaces of known solution
        Eigen::VectorXcd u_t_dir_N = cont_space.Interpolate_helmholtz(u_t_dir,mesh);
        Eigen::VectorXcd u_t_neu_N = cont_space.Interpolate_helmholtz(u_t_neu,mesh);
        Eigen::VectorXcd u_t_N(2*numpanels[i]);
        u_t_N << u_t_dir_N, u_t_neu_N;

        // write difference to computed solution in L^2 norm to file
        filename.open(argv[7], std::ios_base::app);
        filename << mesh.getPanels()[0]->length() << " " << sqrt(abs((sol-u_t_N).dot(M*(sol-u_t_N)))) << std::endl;
        filename.close();
		#ifdef CMDL
        std::cout << "#######################################################" << std::endl;
        std::cout << "Computed Cauchy data on " << numpanels[i] << " panels." << std::endl;
        std::cout << "Residual error of FEM-space interpolation coefficients:" << std::endl;
		std::cout << sqrt(abs((sol-u_t_N).dot(M*(sol-u_t_N)))) << std::endl;
        std::cout << "#######################################################" << std::endl;
        std::cout << std::endl;
		#endif
    }
    return 0;
}
