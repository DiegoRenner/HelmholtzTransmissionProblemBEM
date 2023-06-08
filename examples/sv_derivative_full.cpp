/**
 * \file sv_derivative_full.cpp
 * \brief This target builds a script that computes the singular values and
 * their first two derivatives of the Galerkin BEM
 * approximated BIO for the second-kind direct BIEs of the Helmholtz
 * transmission problem.
 * Minima in the smallest singular value are determined as well
 * by use of the Newton-Raphson method.
 * The scatterer is set to be a circle.
 * The results are written to file.
 * The script can be run as follows:
 *
 * <tt>
 * /path/to/sv_derivative_full \<radius of circle\> \<refraction inside\>
 *      \<refraction outside\> \<initial wavenumber\>
 *      \<\#panels\> \<order of quadrature rule\> \<outputfile\>.
 * </tt>
 *
 * The resulting file will contain the value of \f$k\f$ in the first column.
 * Then the singular values and their first two derivatives at 
 * \f$k\f$ will be listed from smallest to largest in the columns.
 * The singular values and their derivatives occupy three neighboring columns.
 * The final three columns will contain the value of the root,
 * the value of the first derivative at the root and the 
 * number of iterations taken to find the root in the interval 
 * between the current and the next evaluation point.
 * If no root was found these three columns will contain \f$\verb|NAN|\f$.
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

// define shorthand for complex data type and imaginary unit
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

// tolerance when verifying/finding root
double epsilon = 1e-3;//numeric_limits<double>::epsilon();

int main(int argc, char** argv) {

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[2]);
    double c_o = atof(argv[3]);
    complex_t k_0 = atof(argv[4]);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = 250;
    unsigned n_points_y = 1;
    unsigned n_runs_N = 1;
    unsigned numpanels = atoi(argv[5]);
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries
    unsigned order = atoi(argv[6]);
    // define which singular value to evaluate for root finding
    unsigned m = 0;

    // clear existing file
    std::ofstream file_out;
    file_out.open(argv[7], std::ofstream::out | std::ofstream::trunc);
    file_out.close();

	#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Computing singular values and derivatives of BIO." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
	#endif

    // Initialize operators and mesh
    SolutionsOperator so(mesh, order);
    Eigen::MatrixXcd T_next = so.gen_sol_op(k_0, c_o, c_i);
    Eigen::MatrixXcd T_der_next = so.gen_sol_op_1st_der(k_0, c_o, c_i);
    Eigen::MatrixXcd T_der2_next = so.gen_sol_op_2nd_der(k_0, c_o, c_i);
    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            Eigen::MatrixXd res(2*numpanels,3);
            // define wavenumber for current loop
            complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);
            complex_t k_temp_next = (k_0+(j+1)*h_x+ii*double(k)*h_y);

            // compute new operators and assign old ones
            Eigen::MatrixXcd T = T_next;
            Eigen::MatrixXcd T_der = T_der_next;
            Eigen::MatrixXcd T_der2 = T_der2_next;
            T_next = so.gen_sol_op(k_temp_next, c_o, c_i);
            T_der_next = so.gen_sol_op_1st_der(k_temp_next, c_o, c_i);
            T_der2_next = so.gen_sol_op_2nd_der(k_temp_next, c_o, c_i);

            // set singular values and derivatives to compute, all
            unsigned count = T.cols();
            double list[count];
            for (unsigned i = 0; i < count; i++){
                list[i] = i;
            }

            // compute singular value, derivative and 2nd derivative
            res = direct::sv_2nd_der(T,T_der,T_der2,list,count);

            // define functions for computing roots
            auto sv_eval_der = [&] (double k_in) {
                Eigen::MatrixXcd T_in;
                Eigen::MatrixXcd T_der_in;
                if (k_in == k_temp.real()){
                    T_in = T;
                    T_der_in = T_der;
                } else if (k_in == k_temp_next.real()){
                    T_in = T_next;
                    T_der_in = T_der_next;
                } else {
                    T_in = so.gen_sol_op(k_in, c_o, c_i);
                    T_der_in = so.gen_sol_op_1st_der(k_in, c_o, c_i);
                }
                return direct::sv_1st_der(T_in, T_der_in, list, count)(m,1);
            };
            auto sv_eval_both = [&] (double k_in) {
                Eigen::MatrixXcd T_in;
                Eigen::MatrixXcd T_der_in;
                Eigen::MatrixXcd T_der2_in;
                if (k_in == k_temp.real()){
                    T_in = T;
                    T_der_in = T_der;
                    T_der2_in = T_der2;
                } else if (k_in == k_temp_next.real()){
                    T_in = T_next;
                    T_der_in = T_der_next;
                    T_der2_in = T_der2_next;
                } else {
                    T_in = so.gen_sol_op(k_in, c_o, c_i);
                    T_der_in = so.gen_sol_op_1st_der(k_in, c_o, c_i);
                    T_der2_in = so.gen_sol_op_2nd_der(k_in, c_o, c_i);
                }
                Eigen::MatrixXd res = direct::sv_2nd_der(T_in, T_der_in, T_der2_in, list, count).block(m,1,1,2);
                return res;
            };

			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			std::cout << "Singular values and derivatives computed." << std::endl;
			std::cout << "The values for the smallest singular values are:" << std::endl;
			std::cout << res(0, 0) << " " << res(0, 1) << " " << res(0, 2) << std::endl;
			#endif

            // compute minima of smallest singular value
            bool root_found = false;
            unsigned num_iter;
            double root = rtsafe(sv_eval_der, sv_eval_both,k_temp.real(), k_temp.real()+h_x,epsilon,root_found,num_iter);
			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			std::cout << std::endl;
			#endif
            // write results to file for plotting later on
            file_out.open(argv[7], std::ios_base::app);
            file_out << k_temp.real() << " ";
            for (unsigned i = 0; i < count; i++){
                file_out << res(i, 0) << " " << res(i, 1) << " " << res(i, 2) << " ";
            }
            // check that found root is not a crossing
            if (root_found) {
                double val_at_root = sv_eval_der(root);
                if (abs(val_at_root) > epsilon) {
                    file_out << NAN << " " << NAN << " " << NAN << std::endl;
                } else {
                    file_out << root << " " << val_at_root << " " << num_iter << std::endl;
                }
            } else {
                file_out << NAN << " " << NAN << " " << NAN << std::endl;
            }
            file_out.close();
        }
    }
    return 0;
}

