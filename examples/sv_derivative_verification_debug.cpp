
#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "singular_values.hpp"
#include "iomanip"
#include "find_roots.hpp"
#include "single_layer.hpp"
#include "single_layer_der.hpp"
#include "single_layer_der2.hpp"
#include "double_layer.hpp"
#include "double_layer_der.hpp"
#include "double_layer_der2.hpp"
#include "hypersingular.hpp"
#include "hypersingular_der.hpp"
#include "hypersingular_der2.hpp"
#include "discontinuous_space.hpp"
#include "continuous_space.hpp"

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

    // clear existing file
    std::ofstream filename;
    filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/sv_derivative_verification_debug.dat", std::ofstream::out | std::ofstream::trunc);
    filename.close();
    DiscontinuousSpace<0> discont_space;
    ContinuousSpace<1> cont_space;

    // loop over mesh size and wavenumbers
    for (unsigned i = 0; i < n_runs_N; i++) {
        // compute mesh for numpanels
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        for (unsigned j = 0; j < n_points_x; j++) {
            for (unsigned k = 0; k < n_points_y; k++) {
                Eigen::MatrixXd res(2*numpanels[i],3);
                // define wavenumber for current loop
                complex_t k_temp = (k_0+j*h_x+ii*double(k)*h_y);

                //Eigen::MatrixXcd T = single_layer_helmholtz::GalerkinMatrix(mesh,discont_space,order,k_temp,c_i);
                //Eigen::MatrixXcd T_der = single_layer_helmholtz_der::GalerkinMatrix(mesh,discont_space,order,k_temp,c_i);
                //Eigen::MatrixXcd T_der2 = single_layer_helmholtz_der2::GalerkinMatrix(mesh,discont_space,order,k_temp,c_i);
                Eigen::MatrixXcd T = double_layer_helmholtz::GalerkinMatrix(mesh,cont_space,discont_space,order,k_temp,c_i);
                Eigen::MatrixXcd T_der = double_layer_helmholtz_der::GalerkinMatrix(mesh,cont_space,discont_space,order,k_temp,c_i);
                Eigen::MatrixXcd T_der2 = double_layer_helmholtz_der2::GalerkinMatrix(mesh,cont_space,discont_space,order,k_temp,c_i);
                //Eigen::MatrixXcd T = hypersingular_helmholtz::GalerkinMatrix(mesh,cont_space,order,k_temp,c_i);
                //Eigen::MatrixXcd T_der = hypersingular_helmholtz_der::GalerkinMatrix(mesh,cont_space,order,k_temp,c_i);
                //Eigen::MatrixXcd T_der2 = hypersingular_helmholtz_der2::GalerkinMatrix(mesh,cont_space,order,k_temp,c_i);

                unsigned count = 1;
                double list[count];
                for (unsigned i = 0; i < count; i++){
                   list[i] = i;
                }

                // compute singular value, derivative and 2nd derivative
                res = sv_2nd_der(T,T_der,T_der2,list,count);

                // define functions for computing derivatives by extrapolation
                auto sv_eval = [&] (double k_in) {
                    //Eigen::MatrixXcd T_in = single_layer_helmholtz::GalerkinMatrix(mesh,discont_space,order,k_in,c_i);
                    Eigen::MatrixXcd T_in = double_layer_helmholtz::GalerkinMatrix(mesh,cont_space,discont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_in = hypersingular_helmholtz::GalerkinMatrix(mesh,cont_space,order,k_in,c_i);
                    return sv(T_in, list, count)(0);
                };
                auto sv_eval_der = [&] (double k_in) {
                    //Eigen::MatrixXcd T_in = single_layer_helmholtz::GalerkinMatrix(mesh,discont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_der_in = single_layer_helmholtz_der::GalerkinMatrix(mesh,discont_space,order,k_in,c_i);
                    Eigen::MatrixXcd T_in = double_layer_helmholtz::GalerkinMatrix(mesh,cont_space,discont_space,order,k_in,c_i);
                    Eigen::MatrixXcd T_der_in = double_layer_helmholtz_der::GalerkinMatrix(mesh,cont_space,discont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_in = hypersingular_helmholtz::GalerkinMatrix(mesh,cont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_der_in = hypersingular_helmholtz_der::GalerkinMatrix(mesh,cont_space,order,k_in,c_i);
                    return sv_1st_der(T_in, T_der_in, list, count)(m,1);
                };
                auto sv_eval_both = [&] (double k_in) {
                    //Eigen::MatrixXcd T_in = single_layer_helmholtz::GalerkinMatrix(mesh,discont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_der_in = single_layer_helmholtz_der::GalerkinMatrix(mesh,discont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_der2_in = single_layer_helmholtz_der2::GalerkinMatrix(mesh,discont_space,order,k_in,c_i);
                    Eigen::MatrixXcd T_in = double_layer_helmholtz::GalerkinMatrix(mesh,cont_space,discont_space,order,k_in,c_i);
                    Eigen::MatrixXcd T_der_in = double_layer_helmholtz_der::GalerkinMatrix(mesh,cont_space,discont_space,order,k_in,c_i);
                    Eigen::MatrixXcd T_der2_in = double_layer_helmholtz_der2::GalerkinMatrix(mesh,cont_space,discont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_in = hypersingular_helmholtz::GalerkinMatrix(mesh,cont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_der_in = hypersingular_helmholtz_der::GalerkinMatrix(mesh,cont_space,order,k_in,c_i);
                    //Eigen::MatrixXcd T_der2_in = hypersingular_helmholtz_der2::GalerkinMatrix(mesh,cont_space,order,k_in,c_i);
                    return sv_2nd_der(T_in, T_der_in, T_der2_in, list, count).block(m,1,1,2);
                };

                // compute derivatives by extrapolation
                double sv_ext_der1 =  der_by_ext(sv_eval,k_temp.real(),epsilon,epsilon,h_x*epsilon);
                double sv_ext_der2 =  der_by_ext(sv_eval_der,k_temp.real(),epsilon,epsilon,h_x*epsilon);
                //bool root_found = false;
                //double root =  rtsafe(sv_eval_both,k_temp.real(), k_temp.real()+h_x,epsilon,root_found);
                //std::cout << root << std::endl;
                // write results to file for plotting later on
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/sv_derivative_verification_debug.dat", std::ios_base::app);
                filename << k_temp.real() << " "<< res(m,0) << " " << res(m,1)  << " " << sv_ext_der1 << " " << res(m,2) << " " << sv_ext_der2 << std::endl;
                //double val_at_root = sv_eval_der(root);
                //if ( abs(val_at_root) > epsilon) {
                //    filename <<  " " << NAN << " " << NAN << std::endl;
                //} else {
                //    filename <<  " " << root << " " << val_at_root << std::endl;
                //}
                filename.close();
                std::cout << "**********************" << std::endl;


            }
        }
    }
    return 0;
}

