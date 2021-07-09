/**
 * \file roots_newton_circle.cpp
 * \brief This target builds a script that computes minimas in the smallest singular value of the
 * Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
 * transmission problem using the Newton-Raphson method.
 * The scatterer is set to be a circle.
 * The results are written to disk.
 * The script can be run as follows:
 *
 * <tt>
 *  /path/to/roots_newton_circle \<radius of circle\> \<refraction inside\>
 *     \<refraction outside\> \<initial wavenumber\> \<\#grid points for root search\>
 *     \<\#panels\> \<order of quadrature rule\> \<outputfile\>.
 * </tt>
 *
 * The resulting file will contain the left boundary of the 
 * interval used to compute the root in the first column. 
 * Then in the next three columns will be the point, 
 * the function value and the derivative at which the root was found.
 * The last column will contain the number of iterations used to find the root.
 * If no root was found the last four columns will be set to \f$\verb|NAN|\f$.
 * The singular vaues and their derivatives are computed using the direct
 * Eigen algoithm.
 * The user will be updated through the command line about the
 * progress of the algorithm.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#include <complex>
#include <iostream>
#include <fstream>
#include <chrono>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"

// define shorthand for time benchmarking tools, complex data type and immaginary unit
using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

// tolerance when verifying root
double epsilon_ver = 1e-3;
// tolerance when finding root
double epsilon_fin = 1e-6;

int main(int argc, char** argv){

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[2]);
    double c_o = atof(argv[3]);
    complex_t k_0 = atof(argv[4]);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = atoi(argv[5]);
    unsigned n_points_y = 1;
    unsigned numpanels;
    numpanels = atoi(argv[6]);
    double h_x = 100.0/n_points_x;
    double h_y = 100.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = atoi(argv[7]);
    unsigned m = 0;

    // generate output filename with set parameters
    std::string base_order = "../data/file_roots_newton_circle_direct_";
    std::string suffix = ".dat";
    std::string divider = "_";
    std::string file_minimas = base_order.append(argv[2])
                               + suffix;
    // clear existing file
    std::ofstream file_out;
    file_out.open(file_minimas, std::ofstream::out | std::ofstream::trunc);
    file_out.close();

    // Inform user of started computation.
	#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Finding resonances using Newton's method." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
	#endif
    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            auto duration_ops = milliseconds ::zero();
            auto duration = milliseconds::zero();

            // define wavenumber for current loop
            complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

            // set which singular values to evaluate, smallest only
            unsigned count = 1;
            double list[count];
            for (unsigned i = 0; i < count; i++){
                list[i] = i;
            }

            // define functions that return singular value and it's derivative
            auto sv_eval = [&] (double k_in) {
                auto start = high_resolution_clock::now();
                Eigen::MatrixXcd T_in;
                T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                auto end = high_resolution_clock::now();
                duration_ops += duration_cast<milliseconds>(end-start);
                return direct::sv(T_in, list, count)(m);
            };
            auto sv_eval_both = [&] (double k_in) {
                auto start = high_resolution_clock::now();
                Eigen::MatrixXcd T_in;
                Eigen::MatrixXcd T_der_in;
                Eigen::MatrixXcd T_der2_in;
                T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
                T_der2_in = gen_sol_op_2nd_der(mesh, order, k_in , c_o, c_i);
                Eigen::MatrixXd res = direct::sv_2nd_der(T_in, T_der_in, T_der2_in, list, count).block(m,1,1,2);
                auto end = high_resolution_clock::now();
                duration_ops += duration_cast<milliseconds>(end-start);
                return res;
            };
            auto sv_eval_der = [&] (double k_in) {
                auto start = high_resolution_clock::now();
                Eigen::MatrixXcd T_in;
                Eigen::MatrixXcd T_der_in;
                T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
                double res = direct::sv_1st_der(T_in, T_der_in, list, count)(m,1);
                auto end = high_resolution_clock::now();
                duration_ops += duration_cast<milliseconds>(end-start);
                return res;
            };

            // define functions that return singular value and it's derivative
            bool root_found = false;
            unsigned num_iter=0;
            auto start = high_resolution_clock::now();
			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			#endif
            double root =  rtsafe(sv_eval_der,sv_eval_both,k_temp.real(), k_temp.real()+h_x,epsilon_fin,root_found,num_iter);
            auto end = high_resolution_clock::now();
            duration += duration_cast<milliseconds>(end-start);
            // define functions that return singular value and it's derivative
            file_out.open(argv[8], std::ios_base::app);
            file_out << k_temp.real();// << " " << duration.count() << " " << duration_ops.count();
	    
			#ifdef CMDL
	    	// write interval searched to command line
            std::cout << "Interval searched: [" << k_temp.real() 
	    	<< "," << k_temp.real()+h_x << "]" << std::endl;
			#endif

            // check if root was found
            if (root_found) {
                double val_at_root = sv_eval_der(root);
                // check if it's actually a root and not a crossing
                if (abs(val_at_root) < epsilon_ver) {
                    file_out << " " << root << " " << val_at_root << " " << sv_eval(root) << " " << num_iter << std::endl;
				#ifdef CMDL
				// write found roots to command line
				std::cout << "A root was found at: " << root << std::endl;
				std::cout << "The value of the first derivative here is: " << val_at_root << std::endl;
				std::cout << "Number of iterations taken: " << num_iter << std::endl;
				#endif 
                } else {
                    file_out << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
                }
            } else{
                file_out << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
            }
            file_out.close();
			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
            std::cout << std::endl;
			#endif
        }
    }
    return 0;
}

