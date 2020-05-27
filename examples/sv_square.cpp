
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_line.hpp"
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
    using PanelVector = PanelVector;
  // Corner points for the square
  Eigen::RowVectorXd x1(2);
  x1 << 0,0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << eps, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << eps, eps; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, eps; // Point (0,1.5)
  // Parametrized line segments forming the edges of the polygon
  ParametrizedLine line1(x1, x2);
  ParametrizedLine line2(x2, x3);
  ParametrizedLine line3(x3, x4);
  ParametrizedLine line4(x4, x1);
  // Splitting the parametrized lines into panels for a mesh to be used for
  // BEM (Discretization). Here Split is used with input "1" implying that the
  // original edges are used as panels in our mesh.
  PanelVector line1panels = line1.split(numpanels/4);
  PanelVector line2panels = line2.split(numpanels/4);
  PanelVector line3panels = line3.split(numpanels/4);
  PanelVector line4panels = line4.split(numpanels/4);
  PanelVector panels;
  // Storing all the panels in order so that they form a polygon
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  // Construction of a ParametrizedMesh object from the vector of panels
  ParametrizedMesh parametrizedmesh(panels);


    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = atoi(argv[6]);
    unsigned m = 0;

    // clear existing file
    std::ofstream filename;
    filename.open(argv[7], std::ofstream::out | std::ofstream::trunc);
    filename.close();

    // loop over mesh size and wavenumbers
        // compute mesh for numpanels
        ParametrizedMesh mesh(panels);
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

