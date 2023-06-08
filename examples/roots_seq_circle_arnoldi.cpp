/**
 * \file roots_seq_circle_arnoldi.cpp
 * \brief This target builds a script that computes minimas in the smallest singular value of the
 * Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
 * transmission problem using the Algorithm described in Listing 1 of
 * https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2022/2022-38.pdf.
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
#include "parametrized_circular_arc.hpp"
#include "singular_values_arnoldi.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "operators.hpp"

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
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = atoi(argv[7]);
    unsigned m = 0;

    // define accurracy of arnoldi algorithm
    double acc = atof(argv[8]);

    // generate output filename with set parameters
    std::string base_name = "../data/file_roots_seq_circle_arnoldi_";
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
    std::cout << "Finding resonances using sequential algorithm." << std::endl;
    std::cout << "Computing on userdefined problem using circular domain." << std::endl;
    std::cout << std::endl;
	#endif

    SolutionsOperator so(mesh, order);

    auto duration_ops = milliseconds ::zero();
    auto duration = milliseconds::zero();

    // initialize for function calls counter
    unsigned N_fct_calls = 0;

    // set which singular values to evaluate, smallest only
    unsigned count = 1;

    // define functions that return singular value and it's derivative
    auto sv_eval = [&] (double k_in) {
        auto start = high_resolution_clock::now();
        Eigen::MatrixXcd T_in;
        T_in = so.gen_sol_op(k_in, c_o, c_i);
        auto end = high_resolution_clock::now();
        duration_ops += duration_cast<milliseconds>(end-start);
        return arnoldi::sv(T_in, count, acc)(m);
    };
    auto sv_eval_both = [&] (double k_in) {
        N_fct_calls += 1;
        auto start = high_resolution_clock::now();
        Eigen::MatrixXcd T_in = so.gen_sol_op(k_in, c_o, c_i);
        Eigen::MatrixXcd T_der_in = so.gen_sol_op_1st_der(k_in, c_o, c_i);
        Eigen::MatrixXcd T_der2_in = so.gen_sol_op_2nd_der(k_in, c_o, c_i);
        Eigen::MatrixXd res = arnoldi::sv_2nd_der(T_in, T_der_in, T_der2_in, count, acc).block(m,1,1,2);
        auto end = high_resolution_clock::now();
        duration_ops += duration_cast<milliseconds>(end-start);
        return res;
    };
    auto sv_eval_der = [&] (double k_in) {
        auto start = high_resolution_clock::now();
        Eigen::MatrixXcd T_in = so.gen_sol_op(k_in, c_o, c_i);
        Eigen::MatrixXcd T_der_in = so.gen_sol_op_1st_der(k_in, c_o, c_i);
        double res = arnoldi::sv_1st_der(T_in, T_der_in, count, acc)(m,1);
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
    std::vector<std::vector<grid_data>> records;
    std::function<void(std::vector<grid_data>)> recorder = [&records](std::vector<grid_data> entry)->void{records.push_back(entry);};
    std::vector<double> roots = findZeros_seq(sv_eval_both,k_0.real(),k_0.real()+10.0, n_points_x,recorder);
    auto end = high_resolution_clock::now();
    duration += duration_cast<milliseconds>(end-start);


    // define functions that return singular value and it's derivative
    file_out.open(file_minimas, std::ios_base::app);
    file_out << "Initial grid: " << std::endl;
    file_out << *records.begin() << std::endl;
    file_out << std::endl;
    for (auto it = records.begin()+1; it != records.end(); ++it){
        int i = it-records.begin();
        file_out << "Grid after " << i << " iteration(s): " << std::endl;
        file_out << *it;
        file_out << std::endl;
    }
    file_out << "#function calls: "  << N_fct_calls << std::endl;
    file_out << "Roots found: " << roots << std::endl;
    file_out.close();
#ifdef CMDL
    std::cout << "#######################################################" << std::endl;
    std::cout << std::endl;
    #endif
    return 0;
}

