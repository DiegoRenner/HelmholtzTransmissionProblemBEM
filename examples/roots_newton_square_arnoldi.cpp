/**
 * \file roots_newton_square_arnoldi.cpp
 * \brief This target builds a script that computes minimas in the smallest singular value of the
 * Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
 * transmission problem using the Newton-Raphson method.
 * The scatterer is set to be a square.
 * The results are written to disk.
 * The script can be run as follows:
 *
 * <tt>
 *  /path/to/roots_newton_circle \<half side length of square\> \<refraction inside\>
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
 * The singular values and their derivatives are computed using the Arnoldi algorithm.
 * The user will be updated through the command line about the
 * progress of the algorithm
 * if \f$ \verb|-DCMDL| \f$ is set.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#include <complex>
#include <iostream>
#include <fstream>
#include <chrono>
#include "parametrized_line.hpp"
#include "singular_values_arnoldi.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "continuous_space.hpp"

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
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    // compute mesh for numpanels
    using PanelVector = PanelVector;
    // corner points for the square
    Eigen::RowVectorXd x1(2);
    x1 << 0,0; // point (0,0)
    Eigen::RowVectorXd x2(2);
    x2 << eps, 0; // point (1,0)
    Eigen::RowVectorXd x3(2);
    x3 << eps, eps; // point (1,0.5)
    Eigen::RowVectorXd x4(2);
    x4 << 0, eps; // point (0,1.5)
    // parametrized line segments forming the edges of the polygon
    ParametrizedLine line1(x1, x2);
    ParametrizedLine line2(x2, x3);
    ParametrizedLine line3(x3, x4);
    ParametrizedLine line4(x4, x1);
    // splitting the parametrized lines into panels for a mesh to be used for
    // BEM (Discretization).
    PanelVector line1panels = line1.split(numpanels/4);
    PanelVector line2panels = line2.split(numpanels/4);
    PanelVector line3panels = line3.split(numpanels/4);
    PanelVector line4panels = line4.split(numpanels/4);
    PanelVector panels;
    // storing all the panels in order so that they form a polygon
    panels.insert(panels.end(), line1panels.begin(), line1panels.end());
    panels.insert(panels.end(), line2panels.begin(), line2panels.end());
    panels.insert(panels.end(), line3panels.begin(), line3panels.end());
    panels.insert(panels.end(), line4panels.begin(), line4panels.end());
    // construction of a ParametrizedMesh object from the vector of panels
    ParametrizedMesh mesh(panels);

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = atoi(argv[7]);
    unsigned m = 0;

    // define accurracy of arnoldi algorithm
    double acc = atof(argv[8]);

    // generate output filename with set parameters
    std::string base_name = "../data/file_roots_newton_square_arnoldi_";
    std::string suffix = ".dat";
    std::string divider = "_";
    std::string file_minimas = base_name.append(argv[2]).append(divider).append(argv[5])
                                       .append(divider).append(argv[8]) + suffix;
    // clear existing file
    std::ofstream file_out;
    file_out.open(file_minimas, std::ofstream::out | std::ofstream::trunc);
    file_out.close();

    // Inform user of started computation.
	#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Finding resonances using Newton-Raphson method." << std::endl;
    std::cout << "Computing on userdefined problem using square domain." << std::endl;
    std::cout << std::endl;
	#endif

    ContinuousSpace<1> cont_space;
    SolutionsOperator so(mesh, order, cont_space, cont_space);

    // loop over values of wavenumber
    auto tic = high_resolution_clock::now();
    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            auto duration_ops = milliseconds ::zero();
            auto duration = milliseconds::zero();

            // define wavenumber for current loop
            complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

            // set which singular values to evaluate, smallest only
            unsigned count = 1;

            // define functions that return singular value and it's derivative
            auto sv_eval = [&] (double k_in) {
                auto start = high_resolution_clock::now();
                Eigen::MatrixXcd T_in;
                so.gen_sol_op(k_in, c_o, c_i, T_in);
                auto end = high_resolution_clock::now();
                duration_ops += duration_cast<milliseconds>(end-start);
                return arnoldi::sv(T_in, count, acc)(m);
            };
            auto sv_eval_both = [&] (double k_in) {
                auto start = high_resolution_clock::now();
                Eigen::MatrixXcd T_in, T_der_in, T_der2_in;
                so.gen_sol_op_2nd_der(k_in, c_o, c_i, T_in, T_der_in, T_der2_in);
                Eigen::MatrixXd res = arnoldi::sv_2nd_der(T_in, T_der_in, T_der2_in, count, acc).block(m,1,1,2);
                auto end = high_resolution_clock::now();
                duration_ops += duration_cast<milliseconds>(end-start);
                return res;
            };
            auto sv_eval_der = [&] (double k_in) {
                auto start = high_resolution_clock::now();
                Eigen::MatrixXcd T_in, T_der_in;
                so.gen_sol_op_1st_der(k_in, c_o, c_i, T_in, T_der_in);
                double res = arnoldi::sv_1st_der(T_in, T_der_in, count, acc)(m,1);
                auto end = high_resolution_clock::now();
                duration_ops += duration_cast<milliseconds>(end-start);
                return res;
            };

            bool root_found = false;
            unsigned num_iter=0;
            auto start = high_resolution_clock::now();
	    
			// search for root
			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			#endif
            double root =  rtsafe(sv_eval_der,sv_eval_both,k_temp.real(), k_temp.real()+h_x,epsilon_fin,root_found,num_iter);
            auto end = high_resolution_clock::now();
            duration += duration_cast<milliseconds>(end-start);
	    
			#ifdef CMDL
			// write interval searched to command line
            std::cout << "Interval searched: [" << k_temp.real() 
	    	<< "," << k_temp.real()+h_x << "]" << std::endl;
			#endif

            // write result to file
            file_out.open(file_minimas, std::ios_base::app);
            file_out << k_temp.real();// << " " << duration.count() << " " << duration_ops.count();

            // check if root was found
            if (root_found) {
                Eigen::Vector2d val_at_root = sv_eval_both(root);
                // check if it's actually a root and not a crossing
                if (std::abs(val_at_root[0]) < epsilon_fin) {
                    file_out << " " << root << " " << sv_eval(root) << " " << val_at_root.transpose() << " " << num_iter << std::endl;
				#ifdef CMDL
				// write found roots to command line
				std::cout << "A root was found at: " << root << std::endl;
				std::cout << "The value of the first derivative here is: " << val_at_root << std::endl;
				std::cout << "Number of iterations taken: " << num_iter << std::endl;
				#endif
                } else {
                    file_out << " " << NAN << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
                }
            } else{
                file_out << " " << NAN << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
            }
            file_out.close();
			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			std::cout << std::endl;
			#endif
        }
    }
    auto toc = high_resolution_clock::now();
#ifdef CMDL
    std::cout << "Total time: " << duration_cast<seconds>(toc - tic).count() << std::endl;
#endif
    return 0;
}

