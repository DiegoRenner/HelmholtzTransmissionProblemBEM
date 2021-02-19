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
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "singular_values_arnoldi.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
double epsilon = 1e-3;
extern int iter_counter_restart;
extern int iter_counter_matvec;

int main(int argc, char** argv) {

    double acc = atof(argv[1]);
    // define radius of circle refraction index and initial wavenumber
    double eps = 0.25;
    double c_i = 20.0;
    double c_o = 1.0;
    complex_t k_0 = 0.84;

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = 50;
    unsigned numpanels;
    numpanels = 20;
    double h_x = 0.5/n_points_x;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    ParametrizedMesh mesh(curve.split(numpanels));

    // define order of quadrature rule used to compute matrix entries
    unsigned order = 11;

    // clear existing file
    std::ofstream file_out;
    std::string base_vals = "../data/file_vals_";
    std::string base_roots = "../data/file_roots_";
    std::string base_iter = "../data/file_iter_";
    std::string suffix = ".dat";
    std::string file_vals = base_vals.append(argv[1])+suffix;
    std::string file_roots = base_roots.append(argv[1])+suffix;
    std::string file_iter = base_iter.append(argv[1])+suffix;
    // set singular values to be computed, all
    unsigned count = 1;

    for (unsigned j = 0; j < n_points_x; j++) {
        Eigen::MatrixXd res(2*numpanels,3);
        // define wavenumber for current loop
        complex_t k_temp = (k_0+j*h_x);

        // compute solutions operator
        Eigen::MatrixXcd T = gen_sol_op(mesh, order, k_temp, c_o, c_i);


        auto sv_eval = [&] (double k_in) {
            Eigen::MatrixXcd T_in;
            T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
            double res = arnoldi::sv(T_in, count, acc)(0);
            return res;
        };
        auto sv_eval_der = [&] (double k_in) {
            Eigen::MatrixXcd T_in;
            Eigen::MatrixXcd T_der_in;
            T_in = gen_sol_op(mesh, order, k_in, c_o, c_i);
            T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
            double res = arnoldi::sv_1st_der(T_in, T_der_in, count, acc)(0,1);
            return res;
        };

        // search for root
        bool root_found = false;
        unsigned num_iter = 0;
        double root = zbrent(sv_eval_der,k_temp.real(),
                             k_temp.real()+h_x,epsilon,root_found,num_iter);


        // compute singular value
        res = arnoldi::sv(T,count,acc);

        // write singular values to file
        file_out.open(file_vals, std::ios_base::app);
        file_out << k_temp.real() << " ";
        file_out << res.block(0, 0, count, 1).transpose() << std::endl;
        file_out.close();
        // write result to file
        file_out.open(file_roots, std::ios_base::app);
        file_out << k_temp.real();//<< " " << duration.count() << " " << duration_ops.count();


        // check if root was found
        if (root_found) {
            double val_at_root = sv_eval_der(root);
            // check if it's actually a root and not a crossing
            if (abs(val_at_root) < epsilon) {
                file_out << " " << root << " " << val_at_root << " " << sv_eval(root) << " " << num_iter << std::endl;
                // write found roots to command line
            } else {
                file_out << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
            }
        } else{
            file_out << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
        }
        file_out.close();

        file_out.open(file_iter, std::ios_base::app);
        file_out << k_temp.real() << " " << iter_counter_restart << " "
                    << iter_counter_matvec << std::endl;
        file_out.close();
        iter_counter_matvec = 0;
        iter_counter_restart = 0;
    }
    return 0;
}
