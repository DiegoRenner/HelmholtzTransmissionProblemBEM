/**
 * \file dirichlet.hpp
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem of the form given in \f$\eqref{eq:dirbvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "mass_matrix.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "double_layer.hpp"
#include "hypersingular.hpp"
#include "hypersingular_cont.hpp"
#include "hypersingular_discont.hpp"
#include "parametrized_mesh.hpp"
#include "single_layer.hpp"
#include <fstream>
#include <iostream>
#include <single_layer_smooth.hpp>
#include <single_layer_correction.hpp>

namespace parametricbem2d {
    typedef std::complex<double> complex_t;
    complex_t ii = complex_t(0, 1);
    namespace bvp {
        namespace direct_first_kind {
            Eigen::VectorXcd solve_neumann(const ParametrizedMesh &mesh,
                                             std::function<complex_t(double, double)> u_dir,
                                             std::function<complex_t(double, double)> u_neu,
                                             unsigned int order,
                                             const double k) {
                unsigned int numpanels = mesh.getNumPanels();
                DiscontinuousSpace<0> discont_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                Eigen::VectorXcd u_dir_N = cont_space.Interpolate_helmholtz(u_dir, mesh);
                Eigen::VectorXcd u_neu_N = discont_space.Interpolate_helmholtz(u_neu, mesh);
                Eigen::VectorXcd u_N(2*numpanels);
                u_N << u_dir_N, u_neu_N;
                // Computing V matrix
                std::cout << "starting computation of operators with " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd M =
                        mass_matrix::GalerkinMatrix(mesh,discont_space,cont_space,order);
                std::cout << "mass matrix computed" << std::endl;
                Eigen::MatrixXcd K =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k);
                std::cout << "hypersingular helmholtz computed" << std::endl;
                // Build rhs for solving
                Eigen::VectorXcd rhs = (0.5*M+K.transpose())*u_neu_N;
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(-W);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_problem_L2norm.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_dir_N).dot((sol - u_dir_N)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_problem_Mnorm.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_dir_N).dot(M*(sol - u_dir_N)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_problem_Wnorm.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_dir_N).dot(W*(sol - u_dir_N)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                std::cout << "****************************" << std::endl;
                return sol;
            }

            Eigen::VectorXcd solve_dirichlet(const ParametrizedMesh &mesh,
                                           std::function<complex_t(double, double)> u_dir,
                                           std::function<complex_t(double, double)> u_neu,
                                           unsigned order,
                                           const double k) {
                unsigned int numpanels = mesh.getNumPanels();
                DiscontinuousSpace<0> discont_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                Eigen::VectorXcd u_dir_N = cont_space.Interpolate_helmholtz(u_dir, mesh);
                Eigen::VectorXcd u_neu_N = discont_space.Interpolate_helmholtz(u_neu, mesh);
                Eigen::VectorXcd u_N(2*numpanels);
                u_N << u_dir_N, u_neu_N;
                // Computing V matrix
                std::cout << "starting computation of operators with " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd M =
                        mass_matrix::GalerkinMatrix(mesh,discont_space,cont_space,order);
                std::cout << "mass matrix computed" << std::endl;
                Eigen::MatrixXcd K =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd V_smooth =
                        single_layer_smooth_helmholtz::GalerkinMatrix(mesh, discont_space, order, k);
                Eigen::MatrixXcd V_correction =
                        single_layer_correction::GalerkinMatrix(mesh, discont_space, order, k);
                Eigen::MatrixXcd V =
                        V_smooth+V_correction;
                std::cout << "single layer helmholtz computed" << std::endl;
                // Build rhs for solving
                Eigen::VectorXcd rhs = (0.5*M-K)*u_dir_N;
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(-V);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_problem_L2norm.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_neu_N).dot((sol - u_neu_N)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_problem_Mnorm.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_neu_N).dot(M*(sol - u_neu_N)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_problem_Vnorm.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_neu_N).dot(V*(sol - u_neu_N)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                std::cout << "****************************" << std::endl;
                return sol;
            }
        } // namespace direct_first_kind
    } // namespace bvp
    namespace tsp {
        namespace direct_second_kind {
            Eigen::VectorXcd solve_debug_3(const ParametrizedMesh &mesh,
                                   std::function<complex_t(double, double)> u_inc_dir,
                                   std::function<complex_t(double, double)> u_inc_neu,
                                   std::function<complex_t(double, double)> sol_dir,
                                   std::function<complex_t(double, double)> sol_neu,
                                   unsigned order,
                                   const double k_o,
                                   const double k_i) {
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd K_o =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_o);
                std::cout << "double layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd K_i =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_i);
                std::cout << "double layer helmholtz computed" << std::endl;
                //std::ifstream fp_data;
                //double real, imag;
                //char sign;
                //int i = 0;
                //std::string path;
                //path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_o_" + std::to_string(numpanels) + ".dat";
                //Eigen::MatrixXcd K_o(numpanels,numpanels);
                //fp_data.open(path);
                //sign='+';
                //while(fp_data >> real >> imag) {
                //    K_o(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                //    i++;
                //    fp_data >> sign >> sign;
                //}
                //fp_data.close();
                //i = 0;
                //Eigen::MatrixXcd K_i(numpanels,numpanels);
                //path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_i_" + std::to_string(numpanels) + ".dat";
                //fp_data.open(path);
                //sign='+';
                //while(fp_data >> real >> imag) {
                //    K_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real,imag);
                //    i++;
                //    fp_data >> sign >> sign;
                //}
                //fp_data.close();
                Eigen::MatrixXcd W_i =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
                std::cout << "hypersingular helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_o =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_o);
                std::cout << "hypersingular helmholtz rhs computed" << std::endl;
                //i = 0;
                //Eigen::MatrixXcd W_i(numpanels,numpanels); // =
                //path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/raw_data/hypersingular_i_" + std::to_string(numpanels) + ".dat";
                //fp_data.open(path);
                //sign='+';
                //while(fp_data >> real >> imag) {
                //    W_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                //    i++;
                //    fp_data >> sign >> sign;
                //}
                //fp_data.close();
                //std::cout << "hypersingular helmholtz computed" << std::endl;
                //i = 0;
                //Eigen::MatrixXcd W_o(numpanels,numpanels);
                //path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/raw_data/hypersingular_o_" + std::to_string(numpanels) + ".dat";
                //fp_data.open(path);
                //sign='+';
                //while(fp_data >> real >> imag) {
                //    W_o(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                //    i++;
                //    fp_data >> sign >> sign;
                //}
                //fp_data.close();
                Eigen::MatrixXcd V_o =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_o);
                std::cout << "single layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_i =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_i);
                std::cout << "single layer helmholtz computed" << std::endl;
                double h = mesh.getPanels()[0]->length()/6.;
                Eigen::VectorXcd H = Eigen::VectorXcd::Ones(numpanels)*h;
                Eigen::MatrixXcd M_cont= (4*H).asDiagonal();
                M_cont.diagonal(-1) = H.segment(0,numpanels-1);
                M_cont.diagonal(+1) = H.segment(0,numpanels-1);
                Eigen::MatrixXcd M_discont = 6*(H).asDiagonal();
                std::cout << "mass matrix computed" << std::endl;

                Eigen::MatrixXcd A(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i)+M_cont;
                A.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o-V_i);
                A.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o-W_i;
                A.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        (K_o-K_i).transpose()+M_discont;
                Eigen::MatrixXcd A_o(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A_o.block(0, 0, K_o.rows(), K_o.cols()) = -K_o + 0.5*M_cont;
                A_o.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o;
                A_o.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o;
                A_o.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        K_o.transpose()+0.5*M_discont;
                Eigen::VectorXcd u_inc_dir_N = cont_space.Interpolate_helmholtz(u_inc_dir, mesh);
                Eigen::VectorXcd u_inc_neu_N = discont_space.Interpolate_helmholtz(u_inc_neu, mesh);
                Eigen::VectorXcd sol_dir_N = cont_space.Interpolate_helmholtz(sol_dir, mesh);
                Eigen::VectorXcd sol_neu_N = discont_space.Interpolate_helmholtz(sol_neu, mesh);
                Eigen::VectorXcd u_inc_N(2*numpanels);
                Eigen::VectorXcd u_sol_N(2*numpanels);
                u_inc_N << u_inc_dir_N, u_inc_neu_N;
                u_sol_N << sol_dir_N, sol_neu_N;
                // Build rhs for solving
                Eigen::VectorXcd rhs = (A_o * u_inc_N);
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(A);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::cout << "-----------------"<< std::endl;
                std::cout << sol.segment(0,2*numpanels).transpose() << std::endl;
                std::cout << "**************************"<< std::endl;
                std::cout << u_sol_N.segment(0,2*numpanels).transpose() << std::endl;
                std::cout << "-----------------"<< std::endl;

                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_L2norm_dir.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_sol_N).segment(0, numpanels).dot((sol - u_sol_N).segment(0, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Mnorm_dir.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_sol_N).segment(0, numpanels).dot(M_cont*(sol - u_sol_N).segment(0, numpanels))))  << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Wnorm_dir.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_sol_N).segment(0, numpanels).dot(W_o*(sol - u_sol_N).segment(0, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_L2norm_neu.dat", std::ios_base::app);
                filename << sqrt(abs((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Mnorm_neu.dat", std::ios_base::app);
                filename << sqrt(abs((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot(M_discont.transpose()*(sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels))))  << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Vnorm_neu.dat", std::ios_base::app);
                filename << sqrt(abs((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot(V_o*(sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                std::cout << "****************************" << std::endl;
                return sol;
            }

            Eigen::VectorXcd solve_debug_2(const ParametrizedMesh &mesh,
                                           std::function<complex_t(double, double)> u_inc_dir,
                                           std::function<complex_t(double, double)> u_inc_neu,
                                           std::function<complex_t(double, double)> sol_dir,
                                           std::function<complex_t(double, double)> sol_neu,
                                           unsigned order,
                                           const double k) {
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators" << std::endl;
                std::ifstream fp_data;
                double real, imag;
                char sign;
                int i = 0;
                Eigen::MatrixXcd K_i(numpanels,numpanels);
                std::string path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd V_i(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/single_layer_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag ) {
                    V_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd W_i(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/hypersingular_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    W_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_i_adj(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/adjoint_double_layer_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i_adj(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                Eigen::MatrixXcd A(K_i.rows() + W_i.rows(), K_i.cols() + V_i.cols());
                A.block(0, 0, K_i.rows(), K_i.cols()) =  K_i + 0.5* Eigen::MatrixXcd::Identity(K_i.rows(), K_i.cols());
                A.block(0, K_i.cols(), V_i.rows(), V_i.cols()) = -V_i;
                A.block(K_i.rows(), 0, W_i.rows(), W_i.cols()) = -W_i;
                A.block(K_i.rows(), K_i.cols(), K_i.cols(), K_i.rows()) =
                        -K_i_adj + 0.5 * Eigen::MatrixXcd::Identity(K_i.rows(), K_i.cols());
                Eigen::VectorXcd u_inc_dir_N = cont_space.Interpolate_helmholtz(u_inc_dir, mesh);
                Eigen::VectorXcd u_inc_neu_N = discont_space.Interpolate_helmholtz(u_inc_neu, mesh);
                Eigen::VectorXcd sol_dir_N = cont_space.Interpolate_helmholtz(sol_dir, mesh);
                Eigen::VectorXcd sol_neu_N = discont_space.Interpolate_helmholtz(sol_neu, mesh);
                Eigen::VectorXcd u_inc_N(2*numpanels);
                Eigen::VectorXcd u_sol_N(2*numpanels);
                u_inc_N << u_inc_dir_N, u_inc_neu_N;
                u_sol_N << sol_dir_N, sol_neu_N;
                // Build rhs for solving
                Eigen::VectorXcd rhs = A * u_sol_N;
                std::ofstream filename;
                std::cout << "_________________" << std::endl;
                std::cout << rhs << std::endl;
                std::cout << "*****************" << std::endl;
                std::cout << u_sol_N << std::endl;
                std::cout << "_________________" << std::endl;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_data.dat", std::ios_base::app);
                filename << (2*rhs-u_sol_N).segment(0, numpanels).lpNorm<2>() << " " << sqrt(2 * M_PI / numpanels)
                         << " " << (rhs).segment(0, numpanels).lpNorm<2>() / sqrt(2 * M_PI / numpanels) << " "
                         <<
                         sol_dir_N.cwiseAbs().mean() << std::endl;

                std::ofstream filename1;
                filename1.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_data.dat", std::ios_base::app);
                filename1 << (2*rhs-u_sol_N).segment(numpanels, numpanels).lpNorm<2>() << " "
                          << sqrt(2 * M_PI / numpanels)
                          << " "
                          << (rhs).segment(numpanels, numpanels).lpNorm<2>() / sqrt(2 * M_PI / numpanels)
                          << " " <<
                          sol_neu_N.cwiseAbs().mean() << std::endl;
                return rhs;
            }

            Eigen::VectorXcd solve_debug_1(const ParametrizedMesh &mesh,
                                           std::function<complex_t(double, double)> u_inc_dir,
                                           std::function<complex_t(double, double)> u_inc_neu,
                                           std::function<complex_t(double, double)> sol_dir,
                                           std::function<complex_t(double, double)> sol_neu,
                                           unsigned order,
                                           const double k) {
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators" << std::endl;
                std::ifstream fp_data;
                double real, imag;
                char sign;
                int i = 0;
                Eigen::MatrixXcd K_o(numpanels,numpanels);
                std::string path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_o_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_o(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_i(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_i_adj(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/adjoint_double_layer_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i_adj(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_o_adj(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/adjoint_double_layer_o_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_o_adj(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd V_o(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/single_layer_o_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    V_o(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd V_i(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/single_layer_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag ) {
                    V_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd W_o(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/hypersingular_o_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag ) {
                    W_o(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd W_i(numpanels,numpanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/hypersingular_i_" + std::to_string(numpanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    W_i(i/numpanels,i%numpanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                Eigen::MatrixXcd A(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A.block(0, 0, K_o.rows(), K_o.cols()) = -K_o-K_i;// + Eigen::MatrixXcd::Identity(K_o.rows(), K_o.cols());
                A.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o+V_i;
                A.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o+W_i;
                A.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        (K_o_adj+K_i_adj);// + Eigen::MatrixXcd::Identity(K_o.rows(), K_o.cols());
                Eigen::MatrixXcd A_o(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A_o.block(0, 0, K_o.rows(), K_o.cols()) = -K_o;// + 0.5 * Eigen::MatrixXcd::Identity(K_o.rows(), K_o.cols());
                A_o.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o;
                A_o.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o;
                A_o.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        K_o_adj;// + 0.5 * Eigen::MatrixXcd::Identity(K_o.rows(), K_o.cols());
                Eigen::VectorXcd u_inc_dir_N = cont_space.Interpolate_helmholtz(u_inc_dir, mesh);
                Eigen::VectorXcd u_inc_neu_N = discont_space.Interpolate_helmholtz(u_inc_neu, mesh);
                Eigen::VectorXcd sol_dir_N = cont_space.Interpolate_helmholtz(sol_dir, mesh);
                Eigen::VectorXcd sol_neu_N = discont_space.Interpolate_helmholtz(sol_neu, mesh);
                Eigen::VectorXcd u_inc_N(2*numpanels);
                Eigen::VectorXcd u_sol_N(2*numpanels);
                u_inc_N << u_inc_dir_N, u_inc_neu_N;
                u_sol_N << sol_dir_N, sol_neu_N;
                // Build rhs for solving
                Eigen::VectorXcd rhs = A_o * u_inc_N;
                // Solving for coefficients
                //Eigen::FullPivLU<Eigen::MatrixXd> dec(V);
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(A);
                Eigen::VectorXcd sol = dec.solve(rhs);
                //Eigen::VectorXcd test = rhs-A*u_sol_N;
                std::ofstream filename;
                std::cout << "_________________" << std::endl;
                std::cout << sol << std::endl;
                std::cout << "*****************" << std::endl;
                std::cout << (u_sol_N) << std::endl;
                std::cout << "_________________" << std::endl;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_data.dat", std::ios_base::app);
                filename << (2.0*sol - u_sol_N).segment(0, numpanels).lpNorm<2>() << " " << 2 * M_PI / numpanels
                         << " " << (sol - u_sol_N).segment(0, numpanels).lpNorm<2>() / sqrt(2 * M_PI / numpanels) << " "
                         <<
                         sol_dir_N.cwiseAbs().mean() << std::endl;
                std::ofstream filename1;
                filename1.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_data.dat", std::ios_base::app);
                filename1 << ((2.0*sol.real() - 2.0*ii*sol.imag()) - u_sol_N).segment(numpanels, numpanels).lpNorm<2>() << " "
                          << 2 * M_PI / numpanels
                          << " "
                          << (sol - u_sol_N).segment(numpanels, numpanels).lpNorm<2>() / sqrt(2 * M_PI / numpanels)
                          << " " <<
                          sol_neu_N.cwiseAbs().mean() << std::endl;
                return sol;
            }

            Eigen::VectorXcd solve(const ParametrizedMesh &mesh,
                                   std::function<complex_t(double, double)> u_inc_dir,
                                   std::function<complex_t(double, double)> u_inc_neu,
                                   std::function<complex_t(double, double)> sol_dir,
                                   std::function<complex_t(double, double)> sol_neu,
                                   unsigned order,
                                   const double k_o,
                                   const double k_i) {
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd K_o =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_o);
                std::cout << "double layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd K_i =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_i);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_i =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
                std::cout << "hypersingular helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_o =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k_o);
                std::cout << "hypersingular helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_o =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_o);
                std::cout << "single layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_i =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_i);
                std::cout << "single layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd M_cont =
                        mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
                Eigen::MatrixXcd M_discont =
                        mass_matrix::GalerkinMatrix(mesh,discont_space,discont_space,order);
                std::cout << "mass matrix computed" << std::endl;

                Eigen::MatrixXcd A(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i)+M_cont;
                A.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o-V_i);
                A.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o-W_i;
                A.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        (K_o-K_i).transpose()+M_discont;
                Eigen::MatrixXcd A_o(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A_o.block(0, 0, K_o.rows(), K_o.cols()) = -K_o + 0.5*M_cont;
                A_o.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o;
                A_o.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o;
                A_o.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        K_o.transpose()+0.5*M_discont;
                Eigen::VectorXcd u_inc_dir_N = cont_space.Interpolate_helmholtz(u_inc_dir, mesh);
                Eigen::VectorXcd u_inc_neu_N = discont_space.Interpolate_helmholtz(u_inc_neu, mesh);
                Eigen::VectorXcd sol_dir_N = cont_space.Interpolate_helmholtz(sol_dir, mesh);
                Eigen::VectorXcd sol_neu_N = discont_space.Interpolate_helmholtz(sol_neu, mesh);
                Eigen::VectorXcd u_inc_N(2*numpanels);
                Eigen::VectorXcd u_sol_N(2*numpanels);
                u_inc_N << u_inc_dir_N, u_inc_neu_N;
                u_sol_N << sol_dir_N, sol_neu_N;
                // Build rhs for solving
                Eigen::VectorXcd rhs = (A_o * u_inc_N);
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(A);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::cout << "-----------------"<< std::endl;
                //std::cout << sol.segment(0,2*numpanels).transpose() << std::endl;
                //std::cout << "**************************"<< std::endl;
                //std::cout << u_sol_N.segment(0,2*numpanels).transpose() << std::endl;
                //std::cout << "-----------------"<< std::endl;

                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_L2norm_dir.dat", std::ios_base::app);
                filename << sqrt((sol - u_sol_N).segment(0, numpanels).dot((sol - u_sol_N).segment(0, numpanels))).real() << " "
                         << 2 * M_PI*0.25 / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Mnorm_dir.dat", std::ios_base::app);
                filename << sqrt((sol - u_sol_N).segment(0, numpanels).dot(M_cont*(sol - u_sol_N).segment(0, numpanels))).real()  << " "
                         << 2 * M_PI*0.25 / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Wnorm_dir.dat", std::ios_base::app);
                filename << sqrt((sol - u_sol_N).segment(0, numpanels).dot(W_o*(sol - u_sol_N).segment(0, numpanels))).real() << " "
                         << 2 * M_PI*0.25 / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_L2norm_neu.dat", std::ios_base::app);
                filename << sqrt((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels))).real() << " "
                         << 2 * M_PI*0.25 / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Mnorm_neu.dat", std::ios_base::app);
                filename << sqrt((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot(M_discont.transpose()*(sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels))).real()  << " "
                         << 2 * M_PI*0.25 / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Vnorm_neu.dat", std::ios_base::app);
                filename << sqrt((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot(V_o*(sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels))).real() << " "
                         << 2 * M_PI*0.25 / numpanels << std::endl;
                filename.close();
                //std::cout << "****************************" << std::endl;
                return sol;
            }

            Eigen::VectorXcd solve_smooth(const ParametrizedMesh &mesh,
                                   std::function<complex_t(double, double)> u_inc_dir,
                                   std::function<complex_t(double, double)> u_inc_neu,
                                   std::function<complex_t(double, double)> sol_dir,
                                   std::function<complex_t(double, double)> sol_neu,
                                   unsigned order,
                                   const double k_o,
                                   const double k_i) {
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd K_o =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_o);
                std::cout << "double layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd K_i =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_i);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_i =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
                std::cout << "hypersingular helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_o =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k_o);
                std::cout << "hypersingular helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_o =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_o);
                std::cout << "single layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_i =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_i);
                std::cout << "single layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd M_cont =
                        mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
                Eigen::MatrixXcd M_discont =
                        mass_matrix::GalerkinMatrix(mesh,discont_space,discont_space,order);
                std::cout << "mass matrix computed" << std::endl;

                Eigen::MatrixXcd A(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i)+M_cont;
                A.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o-V_i);
                A.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o-W_i;
                A.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        (K_o-K_i).transpose()+M_discont;
                Eigen::MatrixXcd A_o(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A_o.block(0, 0, K_o.rows(), K_o.cols()) = -K_o + 0.5*M_cont;
                A_o.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o;
                A_o.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o;
                A_o.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        K_o.transpose()+0.5*M_discont;
                Eigen::VectorXcd u_inc_dir_N = cont_space.Interpolate_helmholtz(u_inc_dir, mesh);
                Eigen::VectorXcd u_inc_neu_N = discont_space.Interpolate_helmholtz(u_inc_neu, mesh);
                Eigen::VectorXcd sol_dir_N = cont_space.Interpolate_helmholtz(sol_dir, mesh);
                Eigen::VectorXcd sol_neu_N = discont_space.Interpolate_helmholtz(sol_neu, mesh);
                Eigen::VectorXcd u_inc_N(2*numpanels);
                Eigen::VectorXcd u_sol_N(2*numpanels);
                u_inc_N << u_inc_dir_N, u_inc_neu_N;
                u_sol_N << sol_dir_N, sol_neu_N;
                // Build rhs for solving
                Eigen::VectorXcd rhs = (A_o * u_inc_N);
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(A);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::cout << "-----------------"<< std::endl;
                std::cout << sol.segment(0,2*numpanels).transpose() << std::endl;
                std::cout << "**************************"<< std::endl;
                std::cout << u_sol_N.segment(0,2*numpanels).transpose() << std::endl;
                std::cout << "-----------------"<< std::endl;

                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_L2norm_dir.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_sol_N).segment(0, numpanels).dot((sol - u_sol_N).segment(0, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Mnorm_dir.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_sol_N).segment(0, numpanels).dot(M_cont*(sol - u_sol_N).segment(0, numpanels))))  << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Wnorm_dir.dat", std::ios_base::app);
                filename << sqrt(abs((sol - u_sol_N).segment(0, numpanels).dot(W_o*(sol - u_sol_N).segment(0, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_L2norm_neu.dat", std::ios_base::app);
                filename << sqrt(abs((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Mnorm_neu.dat", std::ios_base::app);
                filename << sqrt(abs((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot(M_discont.transpose()*(sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels))))  << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/transmission_problem_Vnorm_neu.dat", std::ios_base::app);
                filename << sqrt(abs((sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).dot(V_o*(sol.real() + ii*sol.imag() - u_sol_N).segment(numpanels, numpanels)))) << " "
                         << 2 * M_PI / numpanels << std::endl;
                filename.close();
                std::cout << "****************************" << std::endl;
                return sol;
            }

            Eigen::MatrixXcd compute_operator(const ParametrizedMesh &mesh,
                                   unsigned order,
                                   complex_t k_o,
                                   complex_t  k_i) {
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd K_o =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_o);
                std::cout << "double layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd K_i =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k_i);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_i =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k_i);
                std::cout << "hypersingular helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_o =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k_o);
                std::cout << "hypersingular helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_o =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_o);
                std::cout << "single layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_i =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k_i);
                std::cout << "single layer helmholtz computed" << std::endl;

                double h = mesh.getPanels()[0]->length()/6.;
                Eigen::VectorXcd H = Eigen::VectorXcd::Ones(numpanels)*h;
                Eigen::MatrixXcd M_cont= mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
                Eigen::MatrixXcd M_discont = mass_matrix::GalerkinMatrix(mesh,discont_space,discont_space,order);
                Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                M.block(0, 0, K_o.rows(), K_o.cols()) = M_cont;
                M.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) = M_discont;
                std::cout << "mass matrix computed" << std::endl;

                Eigen::MatrixXcd A(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A.block(0, 0, K_o.rows(), K_o.cols()) = -K_o + K_i;
                A.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o-V_i;
                A.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o-W_i;
                A.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        (K_o-K_i).transpose();
                A+=M;
                Eigen::MatrixXcd res(2*numpanels,2*numpanels);
                Eigen::LDLT<Eigen::MatrixXcd> llt(M);
                res.setIdentity();
                res = llt.transpositionsP()*res;
                res = llt.matrixU() * res;
                res = llt.vectorD().cwiseSqrt().asDiagonal() * res;
                Eigen::MatrixXcd L = res.transpose();
                Eigen::HouseholderQR<Eigen::MatrixXcd> qr(L);
                Eigen::MatrixXcd R = qr.matrixQR().triangularView<Eigen::Upper>();
                Eigen::MatrixXcd Q = qr.householderQ();
                Eigen::MatrixXcd T = R.inverse()*Q.transpose();
                return T*A*T.transpose();
            }


        } // namespace direct_second_kind
    } // namespace tsp
} // namespace parametricbem2d
