/**
 * \file sv_derivative_verification_circle.cpp
 * \brief This target builds a script that verifies the derivatives of the singular
 * values of the Galerkin BEM approximated BIO for the
 * second-kind direct BIEs of the Helmholtz transmission problem
 * using extrapolation.
 * The scatterer is set to be a circle.
 * The results are written to file.
 * The script can be run as follows:
 *
 * <tt>
 * /path/to/sv_derivative_verification_circle \<radius of circle\>
 *      \<refraction inside\> \<refraction outside\> \<initial wavenumber\>
 *      \<\#panels\> \<order of quadrature rule\> \<outputfile\>.
 * </tt>
 *
 * The resulting file will contain the value of \f$k\f$ in the first column.
 * The second column will contain the value of the smallest singular value at this \f$k\f$.
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
#include "singular_values.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "continuous_space.hpp"

// define shorthand for complex data type and imaginary unit
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

// tolerance for computing derivative by extrapolation
double epsilon = 1e-6;

int main(int argc, char** argv) {

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[2]);
    double c_o = atof(argv[3]);
    complex_t k_0 = atof(argv[4]);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = 2500;
    unsigned n_points_y = 1;
    unsigned numpanels;
    numpanels = atoi(argv[5]);
    double h_x = 100.0/n_points_x;
    double h_y = 100.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries
    // and which singular value to evaluate for extrapolation
    unsigned order = atoi(argv[6]);
    unsigned m = 0;

    // clear existing file
    std::ofstream filename;
    filename.open(argv[7], std::ofstream::out | std::ofstream::trunc);
    filename.close();

	#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Validating singular value derivatives of BIO." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
	#endif
    // initialize operators
    ContinuousSpace<1> cont_space;
    SolutionsOperator so(mesh, order, cont_space, cont_space);
    Eigen::MatrixXcd T_next, T_der_next, T_der2_next;
    so.gen_sol_op_2nd_der(k_0, c_o, c_i, T_next, T_der_next, T_der2_next);

    // loop over wavenumber
    for (unsigned j = 0; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                Eigen::MatrixXd res(2*numpanels,3);
                // define wavenumber for current loop
                complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);
                complex_t k_temp_next = (k_0+(j+1)*h_x+ii*double(k)*h_y);

                //compute new operators and reassign old ones
                Eigen::MatrixXcd T = T_next;
                Eigen::MatrixXcd T_der = T_der_next;
                Eigen::MatrixXcd T_der2 = T_der2_next;
                so.gen_sol_op_2nd_der(k_temp_next, c_o, c_i, T_next, T_der_next, T_der2_next);

                // set singular values to be computed, smallest only
                unsigned count = 1;
                double list[count];
                for (unsigned i = 0; i < count; i++){
                   list[i] = i;
                }

                // compute singular value, derivative and 2nd derivative
                res = direct::sv_2nd_der(T,T_der,T_der2,list,count);

                // define functions for computing derivatives by extrapolation
                auto sv_eval = [&] (double k_in) {
                    Eigen::MatrixXcd T_in;
                    if (k_in == k_temp.real()) {
                        T_in = T;
                    } else if (k_in == k_temp_next.real()){
                        T_in = T_next;
                    } else {
                       so.gen_sol_op(k_in, c_o, c_i, T_in);
                    }
                    return direct::sv(T_in, list, count)(m);
                };
                auto sv_eval_der = [&] (double k_in) {
                    Eigen::MatrixXcd T_in;
                    Eigen::MatrixXcd T_der_in;
                    if (k_in == k_temp.real()) {
                        T_in = T;
                        T_der_in = T_der;
                    } else if (k_in == k_temp_next.real()){
                        T_in = T_next;
                        T_der_in = T_der_next;
                    } else {
                        so.gen_sol_op_1st_der(k_in, c_o, c_i, T_in, T_der_in);
                    }
                    return direct::sv_1st_der(T_in, T_der_in, list, count)(m,1);
                };

                // compute derivatives by extrapolation
                double sv_ext_der1 =  direct::der_by_ext(sv_eval,k_temp.real(),epsilon,epsilon,h_x*epsilon);
                double sv_ext_der2 =  direct::der_by_ext(sv_eval_der,k_temp.real(),epsilon,epsilon,h_x*epsilon);

                // write results to file
                filename.open(argv[7], std::ios_base::app);
                filename << k_temp.real() << " "<< res(m,0) << " " << res(m,1)  << " " << sv_ext_der1 << " " << res(m,2) << " " << sv_ext_der2 << std::endl;
                filename.close();
				#ifdef CMDL
				std::cout << "#######################################################" << std::endl;
				std::cout << "SV derivatives validated at " << k_temp.real() << "." << std::endl;
				std::cout << "Computed and extrapolated values are:" << std::endl;
				std::cout <<  res(m,1)  << " " << sv_ext_der1 << " " << res(m,2) << " " << sv_ext_der2 << std::endl;
				std::cout << "#######################################################" << std::endl;
				std::cout << std::endl;
				#endif
            }
        }
    return 0;
}

