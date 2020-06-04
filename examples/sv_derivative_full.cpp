
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
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
    unsigned numpanels[n_runs_N];
    numpanels[0] = 100;
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    for (int i=1; i<n_runs_N; i++) {
        numpanels[i] = 2 * numpanels[i - 1];
    }

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = 11;
    unsigned m = 0;

    // clear existing file
    std::ofstream filename;
    filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/sv_derivative_full.dat", std::ofstream::out | std::ofstream::trunc);
    filename.close();

    // loop over mesh size and wavenumbers
    for (unsigned i = 0; i < n_runs_N; i++) {
        // compute mesh for numpanels
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        Eigen::MatrixXcd T_next = gen_sol_op(mesh, order, k_0 , c_o, c_i);
        Eigen::MatrixXcd T_der_next = gen_sol_op_1st_der(mesh, order, k_0 , c_o, c_i);
        Eigen::MatrixXcd T_der2_next = gen_sol_op_2nd_der(mesh, order, k_0 , c_o, c_i);
        for (unsigned j = 0; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                Eigen::MatrixXd res(2*numpanels[i],3);
                // define wavenumber for current loop
                complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);
                complex_t k_temp_next = (k_0+(j+1)*h_x+ii*double(k)*h_y);

                Eigen::MatrixXcd T = T_next;
                Eigen::MatrixXcd T_der = T_der_next;
                Eigen::MatrixXcd T_der2 = T_der2_next;
                T_next = gen_sol_op(mesh, order, k_temp_next , c_o, c_i);
                T_der_next = gen_sol_op_1st_der(mesh, order, k_temp_next , c_o, c_i);
                T_der2_next = gen_sol_op_2nd_der(mesh, order, k_temp_next , c_o, c_i);

                unsigned count = T.cols();
                double list[count];
                for (unsigned i = 0; i < count; i++){
                   list[i] = i;
                }

                // compute singular value, derivative and 2nd derivative
                res = sv_2nd_der(T,T_der,T_der2,list,count);

                // define functions for computing derivatives by extrapolation
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
                        T_in = gen_sol_op(mesh, order, k_in, c_o, c_i);
                        T_der_in = gen_sol_op_1st_der(mesh, order, k_in, c_o, c_i);
                    }
                    return sv_1st_der(T_in, T_der_in, list, count)(m,1);
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
                        T_in = gen_sol_op(mesh, order, k_in, c_o, c_i);
                        T_der_in = gen_sol_op_1st_der(mesh, order, k_in, c_o, c_i);
                        T_der2_in = gen_sol_op_2nd_der(mesh, order, k_in, c_o, c_i);
                    }
                    //??????????????????????????
                    Eigen::MatrixXd res = sv_2nd_der(T_in, T_der_in, T_der2_in, list, count).block(m,1,1,2);
                    return res;
                };

                // compute derivatives by extrapolation
                bool root_found = false;
                double root =  rtsafe(sv_eval_both,k_temp.real(), k_temp.real()+h_x,epsilon,root_found);
                std::cout << root << std::endl;
                std::cout << root_found << std::endl;
                // write results to file for plotting later on
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/sv_derivative_full.dat", std::ios_base::app);
                filename << k_temp.real() << " ";
                for (unsigned i = 0; i < count; i++){
                    filename << res(i,0) << " " << res(i,1)  << " " << res(i,2) << " ";
                }
                if (root_found) {
                    double val_at_root = sv_eval_der(root);
                    if (abs(val_at_root) > epsilon) {
                        filename << NAN << " " << NAN << std::endl;
                    } else {
                        filename << root << " " << val_at_root << std::endl;
                    }
                } else {
                    filename << NAN << " " << NAN << std::endl;
                }
                filename.close();
                std::cout << "**********************" << std::endl;


            }
        }
    }
    return 0;
}

