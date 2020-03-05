#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "solvers.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main() {

    // define radius of circle
    double eps = 0.25;
    double n_i = 23.;
    double n_o = 1.;
    double k_0 = 1.0;
    unsigned n_runs = 1000;
    unsigned numpanels = 100;
    parametricbem2d::ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    unsigned order = 11;
    // Loop over number of panels
    double k_o;
    double k_i;
    for (unsigned i = 0; i < n_runs; i++) {
        k_o = (k_0+i*0.0001)*sqrt(n_o);
        k_i = (k_0+i*0.0001)*sqrt(n_i);
        parametricbem2d::ParametrizedMesh mesh(curve.split(numpanels));
        std::cout << "*****************************" << std::endl;
        Eigen::MatrixXcd A = parametricbem2d::tsp::direct_second_kind::compute_operator(
                mesh, order, k_o, k_i);
        Eigen::BDCSVD<Eigen::MatrixXcd> svd(A);
        std::cout << svd.singularValues().minCoeff() << " " << k_o << std::endl;
        std::ofstream filename;
        filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/SV_analysis.dat", std::ios_base::app);
        filename << svd.singularValues().minCoeff() << " " << k_o << std::endl;
        filename.close();
        std::cout << "*****************************" << std::endl;
    }
    return 0;
}
