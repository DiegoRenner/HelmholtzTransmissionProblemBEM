#include <complex>
#include <iostream>
#include <fstream>
#include <solutions_test.hpp>
#include "parametrized_circular_arc.hpp"
#include "compute_SV_der.hpp"
#include "iomanip"
#include "find_zeros.hpp"

using namespace std;

double epsilon = 1e-6;//numeric_limits<double>::epsilon();
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {

    // define radius of circle
    double eps = 0.25;
    double n_i = 23.;
    double n_o = 1.;
    complex_t k_0 = complex_t(0.0,0.0);
    unsigned n_points_x = 100;
    unsigned n_points_y = 1;
    unsigned n_runs_N = 1;
    double numpanels[n_runs_N];
    numpanels[0] = 50;
    double h_x = 10./n_points_x;
    double h_y = 10./n_points_y;
    std::ofstream filename;
    filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis_roots.dat", std::ofstream::out | std::ofstream::trunc);
    filename.close();

    for (int i=1; i<n_runs_N; i++) {
        numpanels[i] = 2 * numpanels[i - 1];
    }
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
                auto sv_der_eval = [&] (double k) {
                    k_o = (k*sqrt(n_o));
                    k_i = (k*sqrt(n_i));
                    return parametricbem2d::eval_SV_der(mesh,order,sqrt(n_i),k_o,k_i).real();
                };
                double zero =  parametricbem2d::zbrent(sv_der_eval, k_temp, k_temp+h_x, epsilon);
                std::cout << zero << std::endl;
                //std::cout << sv_der_eval(k_temp*sqrt(n_o)) << std::endl;
                std::cout << "In interval ["<< k_temp << "," << k_temp+h_x << "]." << endl;
                std::cout << "**********************" << std::endl;

            }
        }
    }
    return 0;
}

