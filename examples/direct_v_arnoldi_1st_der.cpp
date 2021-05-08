/**
 * \file direct_v_arnoldi.cpp
 * \brief This target builds a script that computes
 * the first derivative of the singular values
 * of the Galerkin BEM approximated BIO for the
 * second-kind direct BIEs of the Helmholtz
 * transmission problem, once using the Arnoldi algorithm
 * and once using s direct solver.
 * The scatterer is set to be a circle.
 * The results are written to file.
 * The script can be run as follows:
 *
 * <tt>
 * /path/to/direct_v_arnoldi_1st_der \<radius of circle\> \<number of SV derivatives to be computed\>
 * \<accurracy of arnoldi algorithm\>.
 * </tt>
 *
 * The script will generate four files:
 * file_vals_eig_\<number of SV derivatives to be computed\>_\<accurracy of arnoldi algorithm\>_1stDer.dat,
 * file_vals_arpp_\<number of SV derivatives to be computed\>_\<accurracy of arnoldi algorithm\>_1stDer.dat,
 * file_timings_\<number of SV derivatives to be computed\>_\<accurracy of arnoldi algorithm\>_1stDer.dat,
 * file_iter_\<number of SV derivatives to be computed\>_\<accurracy of arnoldi algorithm\>_1stDer.dat.
 * These will contain the derivatives computed using the direct solver,
 * the derivatives computed using the Arnoldi algorithm,
 * the time taken by the direct solver and the Arnoldi algorithm,
 * and the number of iterations the Arnoldi algorithm took to converge respectively.
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

// define shorthand for time benchmarking tools, complex data type and imaginary unit
using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

// define external variables in order to count iterations of Arnoldi algorithm
extern int iter_counter_restart;
extern int iter_counter_matvec;

int main(int argc, char** argv) {

    //define number of SVs to compute
    unsigned count = atoi(argv[1]);
    double list[count];
    for (unsigned i = 0; i < count; i++){
        list[i] = i;
    }
    //define accuracy for Arnoldi algorithm
    double acc = atof(argv[2]);
    // define radius of circle, refraction indeces and initial wavenumber
    double eps = 0.25;
    double c_i = 20.0;
    double c_o = 1.0;
    complex_t k_0 = 0.0;


    // define mesh on wavenumber and in space on which to perform comparison
    unsigned n_points_x = 100;
    unsigned n_points_y = 1;
    unsigned numpanels = 50;
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries
    unsigned order = 11;

    // define files for writing out later
    std::ofstream file_out;
    std::string base_vals_eig = "../data/file_vals_eig_";
    std::string base_vals_arpp = "../data/file_vals_arpp_";
    std::string base_timings = "../data/file_timings_";
    std::string base_iter = "../data/file_iter_";
    std::string suffix = "_1stDer.dat";
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

    // iterate over mesh on wavenumber
    for (unsigned j = 0; j < n_points_x; j++) {
        for (unsigned k = 0; k < n_points_y; k++) {
            // initialize storage for results
            Eigen::MatrixXd res_direct(2*numpanels,1);
            Eigen::MatrixXd res_arnoldi(2*numpanels,1);

            // define wavenumber for current loop
            complex_t k_temp = k_0+j*h_x+ii*double(k)*h_y;

            // compute solutions operator and derivative
            auto start_op = high_resolution_clock::now();
            Eigen::MatrixXcd T = gen_sol_op(mesh, order, k_temp, c_o, c_i);
            Eigen::MatrixXcd T_der = gen_sol_op_1st_der(mesh, order, k_temp, c_o, c_i);
            auto end_op = high_resolution_clock::now();
            auto duration_op = duration_cast<milliseconds>(end_op-start_op);

            // compute singular value derivatives using the direct solver given by Eigen
            auto start_direct = high_resolution_clock::now();
            res_direct = direct::sv_1st_der(T,T_der,list,count);
            auto end_direct = high_resolution_clock::now();
            auto duration_direct = duration_cast<milliseconds>(end_direct-start_direct);
            // compute singular value using the Arnoldi algorithm
            auto start_arnoldi = high_resolution_clock::now();
            res_arnoldi = arnoldi::sv_1st_der(T, T_der, count,acc);
            auto end_arnoldi = high_resolution_clock::now();
            auto duration_arnoldi = duration_cast<milliseconds>(end_arnoldi-start_arnoldi);

            // write singular values computed by Eigen to file
            file_out.open(file_vals_eig, std::ios_base::app);
            file_out << k_temp.real() << " ";
            file_out << res_direct.block(0, 1, count, 1).transpose() << std::endl;
            file_out.close();
            // write singular values computed by Arnoldi algorithm to file
            file_out.open(file_vals_arpp, std::ios_base::app);
            file_out << k_temp.real() << " ";
            file_out << res_arnoldi.block(0, 1, count, 1).transpose() << std::endl;
            file_out.close();
            // write timings to file
            file_out.open(file_timings, std::ios_base::app);
            file_out << k_temp.real() << " " << duration_direct.count() << " " << duration_arnoldi.count() << " " << duration_op.count() << std::endl;
            file_out.close();
            // write timings to file
            file_out.open(file_iter, std::ios_base::app);
            file_out << k_temp.real() << " " << iter_counter_restart << " "
                     << iter_counter_matvec << std::endl;
            file_out.close();

            // reset external variables for counting in next iteration
            iter_counter_matvec = 0;
            iter_counter_restart = 0;
        }
    }
    return 0;
}
