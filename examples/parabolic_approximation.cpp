/**
 * \file parabolic_approximation.cpp
 * \brief This target builds a script that tries to find minimas in the smallest sinuglar values
 * of the Galerkin BEM approximated solutions operator for the second-kind direct BIEs of 
 * the Helmholtz Transmission problem.
 * The minimas are searched for using a parabolic approximation
 * based on evaluating the smallest singular values and their first
 * two derivatives.
 * The results are written to disk.
 * The script can be run as follows:
 *
 * <tt>
 * /path/to/library/bin/parabolic_approximation \<outfile\>.
 * </tt>
 *
 * In the file the first column contains the initial point used for the parabolic approximation.
 * The next three columns contain the function value and the first 
 * two derivatives at the initial point that were used to compute the parabolic approximation.
 * The user will be updated through the command line about the
 * progress of the algorithm
 * if \f$ \verb|-DCMDL| \f$ is set.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
double epsilon = 1e-3;

int main(int argc, char** argv) {
    // define radius of circle refraction indexces and initial wavenumber
    double eps = 0.25;
    double c_i = 23.;
    double c_o = 1.;
    complex_t k_0 = complex_t(0.75,0.0);
    // define mesh on which to compute BIOs
    unsigned numpanels = 10;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));
    // define order of quadrature rule used to compute matrix entries
    unsigned order = 11;
    // set which singular value to evaluate
    unsigned count = 5;
    double list[count];
    for (unsigned i = 0; i < count; i++){
        list[i] = i;
    }
    // clear existing file
    std::ofstream filename;
    filename.open(argv[1], std::ofstream::out | std::ofstream::trunc);
    filename.close();
    //set range in which to search for next approximated minima
    double step = 0.25;
    complex_t k_temp = k_0;
    // define functions for evaluating singular values and their derivatives
    auto sv_eval = [&] (double k_in) {
        Eigen::MatrixXcd T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
        return sv(T_in, list, count);
    };
    auto sv_eval_der = [&] (double k_in) {
        Eigen::MatrixXcd T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
        Eigen::MatrixXcd T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
        return sv_1st_der(T_in, T_der_in, list, count).block(0,1,count,1);
    };
    auto sv_eval_der2 = [&] (double k_in) {
        Eigen::MatrixXcd T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
        Eigen::MatrixXcd T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
        Eigen::MatrixXcd T_der2_in = gen_sol_op_2nd_der(mesh, order, k_in , c_o, c_i);
        Eigen::MatrixXd res = sv_2nd_der(T_in, T_der_in, T_der2_in, list, count).block(0,2,count,1);
        return res;
    };
    // Inform user of started computation.
	#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Finding resonances using parabolic approximation." << std::endl;
    std::cout << "Computing on baseline problem using " << numpanels << " panels." << std::endl;
    std::cout << std::endl;
	#endif
    // loop over maximum of allowed iterations
    for (unsigned j = 0; j < 1000; j++) {
        // find next approximation of minima using parabolic approximation
        Eigen::VectorXd res =  parabolic_approximation(sv_eval,sv_eval_der,sv_eval_der2,k_temp.real(),step);
        // write values used for parabolic approximation to file for plotting
        filename.open(argv[1], std::ios_base::app);
        filename << k_temp.real() << " " << res.segment(1,3).transpose() << std::endl;
        filename.close();
        double first_der = sv_eval_der(k_temp.real())(0,0);
        // update user
		#ifdef CMDL
        std::cout << "#######################################################" << std::endl;
        std::cout << "The parabolic approximations were computed at: " << k_temp.real() << std::endl;
        std::cout << "The function value and the derivatives used were:" << std::endl;
        std::cout << res.segment(1,3).transpose() << std::endl;
        std::cout << "The current best approximation for a minima is:" << std::endl;
        std::cout << res(0) << std::endl;
        std::cout << "The value of the first derivative at this point is:" << std::endl;
        std::cout << first_der << std::endl;
		#endif
        // if a minima has been found, take a large step to hop out of current basin
        k_temp = res(0);
        if (abs(first_der) < epsilon){
            k_temp+=2*step;
			#ifdef CMDL
            std::cout << "A resonance has been found, jumping to find next." << std::endl;
			#endif
        }
		#ifdef CMDL
        std::cout << "#######################################################" << std::endl;
		std::cout << std::endl;
		#endif
    }
    return 0;
}

