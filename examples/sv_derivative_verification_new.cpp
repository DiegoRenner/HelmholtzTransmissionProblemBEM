
#include <complex>
#include <iostream>
#include <fstream>
#include <solutions_test.hpp>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "iomanip"
#include "roots.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
double epsilon = 1e-3;//numeric_limits<double>::epsilon();
int main() {

    // define radius of circle refraction index and initial wavenumber
    double eps = 0.25;
    double c_i = 23.;
    double c_o = 1.;
    complex_t k_0 = complex_t(0.0,0.0);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_x = 250;
    unsigned n_points_y = 1;
    unsigned n_runs_N = 1;
    unsigned numpanels[n_runs_N];
    numpanels[0] = 50;
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    for (int i=1; i<n_runs_N; i++) {
        numpanels[i] = 2 * numpanels[i - 1];
    }

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = 11;
    unsigned m = 0;

    // clear existing file
    std::ofstream filename;
    //filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/sv_analysis.dat", std::ofstream::out | std::ofstream::trunc);
    //filename.close();

    // loop over mesh size and wavenumbers
    for (unsigned i = 0; i < n_runs_N; i++) {
        // compute mesh for numpanels
        parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels[i]));
        for (unsigned j = 65; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                Eigen::MatrixXd res(2*numpanels[i],3);
                // define wavenumber for current loop
                complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

                // compute singular value, derivative and 2nd derivative
                res = parametricbem2d::sv_2nd_der(mesh, order, k_temp, c_o, c_i);

                // define functions for computing derivatives by extrapolation
                auto sv_eval = [&] (double k) {
                    return parametricbem2d::sv_2nd_der(mesh, order, k, c_o, c_i)(m,0);
                };
                auto sv_eval_der = [&] (double k) {
                    return parametricbem2d::sv_2nd_der(mesh, order, k, c_o, c_i)(m,1);
                };
                auto sv_eval_both = [&] (double k) {
                    return parametricbem2d::sv_2nd_der(mesh, order, k, c_o, c_i).block(m,1,1,2);
                };

                // compute derivatives by extrapolation
                double sv_ext_der2 =  parametricbem2d::der_by_ext(sv_eval_der,k_temp.real(),epsilon,epsilon,h_x*epsilon);
                double sv_ext_der1 =  parametricbem2d::der_by_ext(sv_eval,k_temp.real(),epsilon,epsilon,h_x*epsilon);
                double root =  parametricbem2d::rtsafe(sv_eval_both,k_temp.real(), k_temp.real()+h_x,epsilon);
                std::cout << root << std::endl;
                // write results to file for plotting later on
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/sv_analysis.dat", std::ios_base::app);
                filename << k_temp.real() << " "<< res(m,0) << " " << res(m,1)  << " " << sv_ext_der1 << " " << res(m,2) << " " << sv_ext_der2;
                double val_at_root = sv_eval_der(root);
                if ( abs(val_at_root) > epsilon) {
                    filename <<  " " << NAN << " " << NAN << std::endl;
                } else {
                    filename <<  " " << root << " " << val_at_root << std::endl;
                }
                filename.close();
                std::cout << "**********************" << std::endl;


            }
        }
    }
    return 0;
}

