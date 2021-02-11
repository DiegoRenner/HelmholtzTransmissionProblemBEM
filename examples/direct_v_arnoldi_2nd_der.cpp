/**
 * \file sv_circle.cpp
 * \brief This target builds a script that computes the singular values
 * of the Galerkin BEM approximated BIO for the
 * second-kind direct BIEs of the Helmholtz
 * transmission problem.
 * The scatterer is set to be a circle.
 * The results are written to file.
 * The script can be run as follows:
 *
 * <tt>
 * /path/to/sv_circle \<radius of circle\> \<refraction inside\>
 *      \<refraction outside\> \<initial wavenumber\>
 *      \<\#panels\> \<order of quadrature rule\> \<outputfile\>.
 * </tt>
 *
 * The resulting file will contain the value of \f$k\f$ in the first column.
 * The rest of the columns contain the singular values from 
 * smallest to largest for this \f$k\f$.
 * The user will be updated through the command line about the
 * progress of the algorithm
 * if \f$ \verb|-DCMDL| \f$ is set.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */
#include <complex>
#include <iostream>
#include <fstream>
#include <iostream>
#include <chrono>
#include "parametrized_circular_arc.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "singular_values.hpp"
#include "singular_values_arnoldi.hpp"

using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
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
    unsigned order = atoi(argv[6]);

    // clear existing file
    std::ofstream file_out;
    file_out.open(argv[7], std::ofstream::out | std::ofstream::trunc);
    file_out.close();
    file_out.open(argv[8], std::ofstream::out | std::ofstream::trunc);
    file_out.close();
    file_out.open(argv[9], std::ofstream::out | std::ofstream::trunc);
    file_out.close();

    // Inform user of started computation.
	#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Computing singular values of BIO." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
	#endif
    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            // initialize storage for results
            Eigen::MatrixXd res_direct(2*numpanels,1);
            Eigen::MatrixXd res_arnoldi(2*numpanels,1);
            // define wavenumber for current loop
            complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

            // compute solutions operator
            auto start_op = high_resolution_clock::now();
            Eigen::MatrixXcd T = gen_sol_op(mesh, order, k_temp, c_o, c_i);
            Eigen::MatrixXcd T_der = gen_sol_op_1st_der(mesh, order, k_temp, c_o, c_i);
            Eigen::MatrixXcd T_der2 = gen_sol_op_2nd_der(mesh, order, k_temp, c_o, c_i);
            auto end_op = high_resolution_clock::now();
            auto duration_op = duration_cast<milliseconds>(end_op-start_op);

            // set singular values to be computed, all
            unsigned count = 5;//T.cols();
            double list[count];
            for (unsigned i = 0; i < count; i++){
                list[i] = i;
            }

            // compute singular value using the direct solver given by Eigen
            auto start_direct = high_resolution_clock::now();
            res_direct = direct::sv_2nd_der(T,T_der,T_der2,list,count);
            auto end_direct = high_resolution_clock::now();
            auto duration_direct = duration_cast<milliseconds>(end_direct-start_direct);
            // compute singular value using the Arnoldi algorithm
            auto start_arnoldi = high_resolution_clock::now();
            res_arnoldi = arnoldi::sv_2nd_der(T,T_der,T_der2,count);
            auto end_arnoldi = high_resolution_clock::now();
            auto duration_arnoldi = duration_cast<milliseconds>(end_arnoldi-start_arnoldi);

            // write singular values computed by Eigen to file
            file_out.open(argv[7], std::ios_base::app);
            file_out << k_temp.real() << " ";
            file_out << res_direct.block(0, 0, count, 1).transpose() << std::endl;
            file_out.close();
            std::cout << res_direct.transpose() << std::endl;
            std::cout << "********************" << std::endl;
            // write singular values computed by Arnoldi algorithm to file
            file_out.open(argv[8], std::ios_base::app);
            file_out << k_temp.real() << " ";
            file_out << res_arnoldi.block(0, 0, count, 1).transpose() << std::endl;
            file_out.close();
            std::cout << res_arnoldi.transpose()<< std::endl;
            std::cout << "********************" << std::endl;
            // write timings to file
            file_out.open(argv[9], std::ios_base::app);
            file_out << k_temp.real() << " " << duration_direct.count() << " " << duration_arnoldi.count() << " " << duration_op.count() << std::endl;
            file_out.close();

			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			std::cout << "Singular values at " << k_temp << " computed." << std::endl;
			std::cout << "Smallest singular value is: "
				<< res.block(0, 0, 1, 1).transpose() << std::endl;
            std::cout << "#######################################################" << std::endl;
			std::cout << std::endl;
			#endif
        }
    }
    return 0;
}