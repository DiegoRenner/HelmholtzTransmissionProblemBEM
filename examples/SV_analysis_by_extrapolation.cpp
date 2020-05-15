#include <complex>
#include <iostream>
#include <fstream>
#include <solutions_test.hpp>
#include "parametrized_circular_arc.hpp"
#include "compute_SV_der.hpp"
#include "iomanip"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {

    // define radius of circle
    double eps = 0.25;
    double n_i = 23.;
    double n_o = 1.;
    complex_t k_0 = complex_t(0.0,0.0);
    unsigned n_points_x = 250;
    unsigned n_points_y = 1;
    unsigned n_runs_N = 1;
    double numpanels[n_runs_N];
    numpanels[0] = 10;
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;

    for (int i=1; i<n_runs_N; i++) {
        numpanels[i] = 2 * numpanels[i - 1];
    }
    std::ofstream filename;
    filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_ext.dat", std::ofstream::out | std::ofstream::trunc);
    filename.close();
    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    unsigned order = 11;
    // Loop over number of panels
    double k_o;
    double k_i;
    for (unsigned i = 0; i < n_runs_N; i++) {
        parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels[i]));
        for (unsigned j = 0; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                double k_temp = (k_0+j*h_x+ii*double(k)*h_y).real();
                auto sv_eval = [&] (double k) {
                    k_o = (k*sqrt(n_o));
                    k_i = (k*sqrt(n_i));
                    return parametricbem2d::eval_SV_der(mesh,order,sqrt(n_i),k_o,k_i).real();
                };
                double der =  parametricbem2d::diffex(sv_eval,k_temp,0.02,0.000001,0.0000000001);
                std::cout << "**********************" << std::endl;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_ext.dat", std::ios_base::app);
                filename << der << " " << sv_eval(k_temp) << " " << k_temp << std::endl;
                filename.close();

            }
        }
    }
    return 0;
}
