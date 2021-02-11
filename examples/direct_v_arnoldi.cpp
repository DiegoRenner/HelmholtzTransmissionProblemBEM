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
#include <chrono>
#include "parametrized_circular_arc.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "singular_values.hpp"
#include "singular_values_arnoldi.hpp"
#include "st_vec_storage.hpp"

using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
extern int iter_counter_restart;
extern int iter_counter_matvec;

int main(int argc, char** argv) {

    double acc = atof(argv[2]);
    // define radius of circle refraction index and initial wavenumber
    double eps = 0.25;
    double c_i = 20.0;
    double c_o = 1.0;
    complex_t k_0 = 0.0;


    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = 100;
    unsigned n_points_y = 1;
    unsigned numpanels;
    numpanels = 50;
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries
    unsigned order = 11;

    // clear existing file
    std::ofstream file_out;

    // Inform user of started computation.
	/*#ifdef CMDL
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "Computing singular values of BIO." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
    #endif*/
    unsigned count = atoi(argv[1]);//T.cols();
    double list[count];
    for (unsigned i = 0; i < count; i++){
        list[i] = i;
    }
    std::string base_vals_eig = "file_vals_eig_";
    std::string base_vals_arpp = "file_vals_arpp_";
    std::string base_timings = "file_timings_";
    std::string base_iter = "file_iter_";
    std::string suffix = ".dat";
    std::string file_vals_eig = base_vals_eig.append(argv[1]);
    std::string file_vals_arpp = base_vals_arpp.append(argv[1]);
    std::string file_timings = base_timings.append(argv[1]);
    std::string file_iter = base_iter.append(argv[1]);
    file_vals_eig = file_vals_eig.append("_");
    file_vals_arpp = file_vals_arpp.append("_");
    file_timings = file_timings.append("_");
    file_iter = file_iter.append("_");
    file_vals_eig = file_vals_eig.append(argv[2])+suffix;
    file_vals_arpp = file_vals_arpp.append(argv[2])+suffix;
    file_timings = file_timings.append(argv[2])+suffix;
    file_iter = file_iter.append(argv[2])+suffix;
    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            // initialize storage for results
            Eigen::MatrixXd res_direct(2*numpanels,1);
            Eigen::MatrixXd res_arnoldi(2*numpanels,1);
            // define wavenumber for current loop
            complex_t k_temp = k_0+j*h_x+ii*double(k)*h_y;

            // compute solutions operator
            auto start_op = high_resolution_clock::now();
            Eigen::MatrixXcd T = gen_sol_op(mesh, order, k_temp, c_o, c_i);
            auto end_op = high_resolution_clock::now();
            auto duration_op = duration_cast<milliseconds>(end_op-start_op);

            // set singular values to be computed, all

            // compute singular value using the direct solver given by Eigen
            auto start_direct = high_resolution_clock::now();
            res_direct = direct::sv(T,list,count);
            auto end_direct = high_resolution_clock::now();
            auto duration_direct = duration_cast<milliseconds>(end_direct-start_direct);
            // compute singular value using the Arnoldi algorithm
            auto start_arnoldi = high_resolution_clock::now();
            res_arnoldi = arnoldi::sv(T, count, acc);
            auto end_arnoldi = high_resolution_clock::now();
            auto duration_arnoldi = duration_cast<milliseconds>(end_arnoldi-start_arnoldi);

            // write singular values computed by Eigen to file
            file_out.open(file_vals_eig, std::ios_base::app);
            file_out << k_temp.real() << " ";
            file_out << res_direct.block(0, 0, count, 1).transpose() << std::endl;
            file_out.close();
            // write singular values computed by Arnoldi algorithm to file
            file_out.open(file_vals_arpp, std::ios_base::app);
            file_out << k_temp.real() << " ";
            file_out << res_arnoldi.block(0, 0, count, 1).transpose() << std::endl;
            file_out.close();
            // write timings to file
            file_out.open(file_timings, std::ios_base::app);
            file_out << k_temp.real() << " " << duration_direct.count() << " " << duration_arnoldi.count() << " " << duration_op.count() << std::endl;
            file_out.close();

            file_out.open(file_iter, std::ios_base::app);
            file_out << k_temp.real() << " " << iter_counter_restart << " "
                     << iter_counter_matvec << std::endl;
            iter_counter_matvec = 0;
            iter_counter_restart = 0;

			#ifdef CMDL
            std::cout << "#######################################################" << std::endl;
			std::cout << "Singular values at " << k_temp << " computed." << std::endl;
			std::cout << "Smallest singular value is: "
				<< res_direct.block(0, 0, 1, 1).transpose() << std::endl;
            std::cout << "#######################################################" << std::endl;
			std::cout << std::endl;
			#endif
        }
    }
    return 0;
}
