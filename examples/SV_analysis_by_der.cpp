#include <complex>
#include <iostream>
#include <fstream>
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
    complex_t k_0 = complex_t(0.1,0.0);
    unsigned n_points_x = 500;
    unsigned n_points_y = 1;
    unsigned n_runs_N = 1;
    double numpanels[n_runs_N];
    numpanels[0] = 50;

    for (int i=1; i<n_runs_N; i++) {
        numpanels[i] = 2 * numpanels[i - 1];
    }
    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    unsigned order = 11;
    // Loop over number of panels
    double k_o = (k_0*sqrt(n_o)).real();
    double k_i = (k_0*sqrt(n_i)).real();
    for (unsigned i = 0; i < n_runs_N; i++) {
        parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels[i]));
        for (unsigned j = 0; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                k_o = ((k_0+j*9.9/n_points_x+ii*double(k)*2.0/double(n_points_y))*sqrt(n_o)).real();
                k_i = ((k_0+j*9.9/n_points_x+ii*double(k)*2.0/double(n_points_y))*sqrt(n_i)).real();
                std::cout << "*****************************" << std::endl;
                complex_t res = parametricbem2d::compute_SV_der(mesh, order, k_o, k_i);
            }
        }
    }
    return 0;
}
