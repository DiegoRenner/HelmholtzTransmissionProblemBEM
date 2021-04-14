/**
 * \file transmission_problem_verification.cpp
 * \brief This target builds a script that computes solutions to
 * the analytically solvable case of the Helmholtz transmission
 * problem where the scatterer is a circle using second-kind direct
 * BIEs and Galerkin BEM.
 * The results are written to file.
 * The script can be run as follows:
 * <tt>
 * /path/to/transmission_problem_verification \<radius of circle\>
 *      \<\#coeffs for series expansion of solution\> \<refraction inside\>
 *      \<refraction outside\> \<wavenumber\>
 *      \<order of quadrature rule\> \<outputfile\>
 * </tt>
 * This output file will contain two columns.
 * The first will contain the current panel size.
 * The second will contain the residual error in the euclidean 
 * norm of the computed FEM-space interpolation coefficients to the 
 * known FEM-space interpolation coefficients for the current number of panels.
 * Then the columns will contain the computed derivative, the 
 * extrapolated derivative, the computed second derivative and the extrapolated second 
 * derivative in this order.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "solvers.hpp"
#include "continuous_space.hpp"
#include "single_layer.hpp"
#include "hypersingular.hpp"
#include "double_layer.hpp"


typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main(int argc, char** argv) {

    // define radius of circle refraction index and initial wavenumber
    double eps = 1;
    double c_o = 1.0;
    std::vector<double> k = {0.5, 1.0, 1.5};

    //set boundary
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);

    // define order of quadrature rule used to compute matrix entries
    unsigned order = 11;


    // set number of panels for which to compute BIOs
    unsigned n_runs = 1;
    int numpanels[n_runs];
    numpanels[0] = 50;
    for (int i=1; i<n_runs; i++){
        numpanels[i] = 2*numpanels[i-1];
    }
    // define incoming wave and resulting wave
    auto par_exp = [&] (double x1, double x2, int n) {
        double t = x2 > 0 ? acos(x1) : 2*M_PI-acos(x1);
        double c = 1/sqrt(2*M_PI);
        return c*exp(ii*t*double(n));
    };
    ContinuousSpace<1> cont_space;

    std::string file_SV_debugging = "../data/file_SV_debugging.txt";
    std::ofstream file_out;
    file_out.open(file_SV_debugging, std::ofstream::out | std::ofstream::trunc);
    file_out.close();
    file_out.open(file_SV_debugging, std::ofstream::app);

    for (int i = 0; i < n_runs; i ++){
        file_out << "#panels = " << numpanels[i] << std::endl;
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        file_out << "Single Layer BIO" << std::endl;
        file_out << "k = [" << k[0] << " " << k[1] << " " << k[2] << "]" <<std::endl;
        for (int j = 0; j < 10; j++) {
            auto par_exp_n = [&](double x1, double x2) {
                int n = j + 1;
                return par_exp(x1, x2, n);
            };
            file_out << "n = " << j << " : ";
            for (int l = 0; l < 3; l++) {
                Eigen::MatrixXcd V =
                        single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k[l], c_o);
                Eigen::VectorXcd par_exp_N = cont_space.Interpolate_helmholtz(par_exp_n,mesh);
                file_out << par_exp_N.conjugate().transpose()*V*par_exp_N << " ";
            }
            file_out << std::endl;
        }
        file_out << std::endl;
        file_out << "Double Layer BIO" << std::endl;
        file_out << "k = [" << k[0] << " " << k[1] << " " << k[2] << "]" <<std::endl;
        for (int j = 0; j < 10; j++) {
            auto par_exp_n = [&](double x1, double x2) {
                int n = j + 1;
                return par_exp(x1, x2, n);
            };
            file_out << "n = " << j << " : ";
            for (int l = 0; l < 3; l++) {
                Eigen::MatrixXcd K =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k[l], c_o);
                Eigen::VectorXcd par_exp_N = cont_space.Interpolate_helmholtz(par_exp_n,mesh);
                file_out << par_exp_N.conjugate().transpose()*K*par_exp_N << " ";
            }
            file_out << std::endl;
        }
        file_out << std::endl;
        file_out << "Hypersingular BIO" << std::endl;
        file_out << "k = [" << k[0] << " " << k[1] << " " << k[2] << "]" <<std::endl;
        for (int j = 0; j < 10; j++) {
            auto par_exp_n = [&](double x1, double x2) {
                int n = j + 1;
                return par_exp(x1, x2, n);
            };
            file_out << "n = " << j << " : ";
            for (int l = 0; l < 3; l++) {
                Eigen::MatrixXcd W =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k[l], c_o);
                Eigen::VectorXcd par_exp_N = cont_space.Interpolate_helmholtz(par_exp_n,mesh);
                file_out << par_exp_N.conjugate().transpose()*W*par_exp_N << " ";
            }
            file_out << std::endl;
        }
        file_out << std::endl;
    }
    file_out.close();

    // generate outputfilename


    return 0;
}
