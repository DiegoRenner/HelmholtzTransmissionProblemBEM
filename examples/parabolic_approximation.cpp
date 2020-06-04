
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "iomanip"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
double epsilon = 1e-3;//numeric_limits<double>::epsilon();
int main() {

    // define radius of circle refraction index and initial wavenumber
    double eps = 0.25;
    double c_i = 23.;
    double c_o = 1.;
    complex_t k_0 = complex_t(0.75,0.0);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = 250;
    unsigned n_points_y = 1;
    unsigned n_runs_N = 1;
    unsigned numpanels[n_runs_N];
    numpanels[0] = 10;
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    for (int i=1; i<n_runs_N; i++) {
        numpanels[i] = 2 * numpanels[i - 1];
    }

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = 11;
    unsigned m = 0;
    std::ofstream filename;
    filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/parabolic_approximation.dat", std::ofstream::out | std::ofstream::trunc);
    filename.close();

    double step = 0.25;
    // clear existing file

    // loop over mesh size and wavenumbers
    for (unsigned i = 0; i < n_runs_N; i++) {
        // compute mesh for numpanels
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        complex_t k_temp = k_0;
        for (unsigned j = 0; j < 1000; j++) {
                // define wavenumber for current loop


                unsigned count = 5;
                double list[count];
                for (unsigned i = 0; i < count; i++){
                   list[i] = i;
                }


                // define functions for computing derivatives by extrapolation
                auto sv_eval = [&] (double k_in) {
                    Eigen::MatrixXcd T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                    return sv(T_in, list, count);
                };
                auto sv_eval_der = [&] (double k_in) {
                    Eigen::MatrixXcd T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                    Eigen::MatrixXcd T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
                    return sv_1st_der(T_in, T_der_in, list, count).block(0,1,count,1);
                };
                auto sv_eval_der2 = [&] (double k_in) {
                    Eigen::MatrixXcd T_in = gen_sol_op(mesh, order, k_in , c_o, c_i);
                    Eigen::MatrixXcd T_der_in = gen_sol_op_1st_der(mesh, order, k_in , c_o, c_i);
                    Eigen::MatrixXcd T_der2_in = gen_sol_op_2nd_der(mesh, order, k_in , c_o, c_i);
                    Eigen::MatrixXd res = sv_2nd_der(T_in, T_der_in, T_der2_in, list, count).block(0,2,count,1);
                    return res;
                };

                // compute derivatives by extrapolation
            Eigen::VectorXd res =  parabolic_approximation(sv_eval,sv_eval_der,sv_eval_der2,k_temp.real(),step);
                    filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/parabolic_approximation.dat", std::ios_base::app);
            filename << k_temp.real() << " " << res.segment(1,3).transpose() << std::endl;
                filename.close();
                k_temp = res(0);
                if (abs(res(2)) < epsilon){
                    k_temp+=2*step;
                }
            std::cout << "**********************" << std::endl;


            }
        }
    return 0;
    }

