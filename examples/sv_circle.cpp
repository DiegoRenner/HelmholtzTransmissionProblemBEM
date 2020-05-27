
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "roots.hpp"
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
    unsigned numpanels;
    numpanels = atoi(argv[5]);
    double h_x = 10.0/n_points_x;
    double h_y = 10.0/n_points_y;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = atoi(argv[6]);
    unsigned m = 0;

    // clear existing file
    std::ofstream filename;
    filename.open(argv[7], std::ofstream::out | std::ofstream::trunc);
    filename.close();

    // loop over mesh size and wavenumbers
        // compute mesh for numpanels
        ParametrizedMesh mesh(curve.split(numpanels));
        for (unsigned j = 0; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                Eigen::MatrixXd res(2*numpanels,3);
                // define wavenumber for current loop
                complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

                Eigen::MatrixXcd T = gen_sol_op(mesh, order, k_temp, c_o, c_i);

                unsigned count = T.cols();
                double list[count];
                for (unsigned i = 0; i < count; i++){
                   list[i] = i;
                }

                // compute singular value, derivative and 2nd derivative
                res = sv(T,list,count);

                // define functions for computing derivatives by extrapolation
                filename.open(argv[7], std::ios_base::app);
                filename << k_temp.real() << " ";
                filename << res.block(0,0,count,1).transpose() << std::endl;
                filename.close();
                std::cout << "**********************" << std::endl;


            }
        }

    return 0;
}

