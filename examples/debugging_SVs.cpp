/**
 * \file debugging_SVs.cpp
 * \brief This target builds a script that computes the Eigenvalues
 * of the BIO for the Helmholtz Transmission Problem.
 * The results are written to file.
 * The script can be run as follows:
 * <tt>
 * /path/to/debugging_SVs
 * </tt>
 * The output file will contain a section for each set mesh resolution
 * and each of those sections will contain one section each for every BIO
 * where all Eigenvalues for different wavenumbers will be listed in columns.
 * The Eigenvalues are computed using the facts stated in Lemma 3.22. [TODO: find reference]
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

// defince shorthand for compley data type and immaginary unit
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

int main(int argc, char** argv) {

    // define radius of circle refraction index and wavenumbers to probe
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
    // define basis function of our Hilbert spaces
    auto par_exp = [&] (double x1, double x2, int n) {
        double t = x2 > 0 ? acos(x1) : 2*M_PI-acos(x1);
        double c = 1/sqrt(2*M_PI);
        return c*exp(ii*t*double(n));
    };

    // initialize FEM space
    ContinuousSpace<1> cont_space;

    // initialize file to wirte to
    std::string file_SV_debugging = "../data/file_SV_debugging.dat";
    std::ofstream file_out;
    // clear file if it already exists
    file_out.open(file_SV_debugging, std::ofstream::out | std::ofstream::trunc);
    file_out.close();
    file_out.open(file_SV_debugging, std::ofstream::app);

    // loop over every mesh resolution
    for (int i = 0; i < n_runs; i ++){
        file_out << "#panels = " << numpanels[i] << std::endl;
        ParametrizedMesh mesh(curve.split(numpanels[i]));

        // compute EVs for each BIO at this resolution
        file_out << "Single Layer BIO" << std::endl;
        file_out << "k = [" << k[0] << " " << k[1] << " " << k[2] << "]" <<std::endl;
        for (int j = 0; j < 10; j++) {
            auto par_exp_n = [&](double x1, double x2) {
                int n = j + 1;
                return par_exp(x1, x2, n);
            };
            file_out << "n = " << j+1 << " : ";
            for (int l = 0; l < 3; l++) {
                Eigen::MatrixXcd V =
                        single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k[l], 0., c_o);
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
            file_out << "n = " << j+1 << " : ";
            for (int l = 0; l < 3; l++) {
                Eigen::MatrixXcd K =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k[l], 0., c_o);
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
            file_out << "n = " << j+1 << " : ";
            for (int l = 0; l < 3; l++) {
                Eigen::MatrixXcd W =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k[l], 0., c_o);
                Eigen::VectorXcd par_exp_N = cont_space.Interpolate_helmholtz(par_exp_n,mesh);
                file_out << par_exp_N.conjugate().transpose()*W*par_exp_N << " ";
            }
            file_out << std::endl;
        }
        file_out << std::endl;
    }
    file_out.close();

    return 0;
}
