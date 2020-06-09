
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include <chrono>

using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
double epsilon = 1e-6;//numeric_limits<double>::epsilon();
int main(int argc, char** argv){

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[2]);
    double c_o = atof(argv[3]);
    complex_t k_0 = atof(argv[4]);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = atoi(argv[5]);
    unsigned n_points_y = 1;
    unsigned n_runs_N = 1;
    unsigned numpanels;
    numpanels = atoi(argv[6]);
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = atoi(argv[7]);
    unsigned m = 0;

    // clear existing file
    std::ofstream filename;
    std::cout << argv[8] << std::endl;
    filename.open(argv[8], std::ofstream::out | std::ofstream::trunc);
    filename.close();

    // loop over mesh size and wavenumbers
        // compute mesh for numpanels
        ParametrizedMesh mesh(curve.split(numpanels));
    auto duration_ops = seconds::zero();
    auto duration = seconds::zero();
    for (unsigned j = 0; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                // define wavenumber for current loop
                complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);


                unsigned count = 1;
                double list[count];
                for (unsigned i = 0; i < count; i++){
                   list[i] = i;
                }

                auto sv_eval = [&] (double k_in) {
                    auto start = high_resolution_clock::now();
                    Eigen::MatrixXcd T_in;
                        T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                    auto end = high_resolution_clock::now();
                    duration_ops += duration_cast<seconds>(end-start);
                    return sv(T_in, list, count)(m);
                };
                auto sv_eval_both = [&] (double k_in) {
                    auto start = high_resolution_clock::now();
                    Eigen::MatrixXcd T_in;
                    Eigen::MatrixXcd T_der_in;
                    Eigen::MatrixXcd T_der2_in;
                        T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                        T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
                        T_der2_in = gen_sol_op_2nd_der(mesh, order, k_in , c_o, c_i);
                    //??????????????????????????
                    Eigen::MatrixXd res = sv_2nd_der(T_in, T_der_in, T_der2_in, list, count).block(m,1,1,2);
                    auto end = high_resolution_clock::now();
                    duration_ops += duration_cast<seconds>(end-start);
                    return res;
                };
                auto sv_eval_der = [&] (double k_in) {
                    auto start = high_resolution_clock::now();
                    Eigen::MatrixXcd T_in;
                    Eigen::MatrixXcd T_der_in;
                        T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                        T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
                    auto end = high_resolution_clock::now();
                    duration_ops += duration_cast<seconds>(end-start);
                    return sv_1st_der(T_in, T_der_in, list, count)(m,1);
                };

                // compute derivatives by extrapolation
                bool root_found = false;
                unsigned num_iter;
                double root =  rtsafe(sv_eval_both,k_temp.real(), k_temp.real()+h_x,epsilon,root_found,num_iter);
                std::cout << root << std::endl;
                filename.open(argv[8], std::ios_base::app);
                filename << k_temp.real();
                if (root_found) {
                    double val_at_root = sv_eval_der(root);
                    if (abs(val_at_root) < epsilon) {
                        filename << " " << root << " " << val_at_root << " " << sv_eval(root) << " "  << num_iter << std::endl;
                    } else {
                        filename << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
                    }
                } else{
                    filename << " " << NAN << " " << NAN << " " << NAN << " " << NAN << std::endl;
                }
                filename.close();
                std::cout << "**********************" << std::endl;


            }
        }
    return 0;
}

