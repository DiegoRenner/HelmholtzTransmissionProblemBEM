/**
 * \file dirichlet.hpp
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem of the form given in \f$\eqref{eq:dirbvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef TRANSMISSIONHPP
#define TRANSMISSIONHPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "adj_double_layer.hpp"
#include "mass_matrix.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "double_layer.hpp"
#include "hypersingular.hpp"
#include "integral_gauss.hpp"
#include "parametrized_mesh.hpp"
#include "single_layer.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <fstream>

namespace parametricbem2d {
    namespace transmission_bvp {
/**
 * This namespace contains the solver using the direct first kind method which
 * has the variational formulation as given in \f$\eqref{eq:aVdir}\f$. The
 * Solver uses the lowest order BEM spaces for computation.
 */

        namespace direct_second_kind {
/**
 * This function is used to solve the Dirichlet boundary value problem given
 * in \f$\eqref{eq:dirbvp}\f$ using the variational formulation given in
 * \f$\eqref{eq:l2dv}\f$ after lifting of the variation formulation given in
 * \f$\eqref{eq:bie2nbpvv}\f$ to the \f$L^{2}(\Gamma)\f$ space. The function
 * outputs a vector of estimated Neumann trace of the solution u.
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param g Dirichlet boundary condition in 2D using a function of the form
 *          double(double,double)
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd type representing the Neumann trace of the
 * solution u
 */
            Eigen::VectorXcd solve_4(const ParametrizedMesh &mesh,
                                     std::function<complex<double>(double, double)> u_inc,
                                     std::function<complex<double>(double, double)> u_inc_T,
                                     std::function<complex<double>(double, double)> sol_D,
                                     std::function<complex<double>(double, double)> sol_N,
                                     unsigned order, const double c_o, const double c_i) {
                typedef complex<double> complex_t;
                complex_t ii = complex_t(0, 1);
                unsigned int numpanels = mesh.getNumPanels();
                DiscontinuousSpace<0> discont_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                Eigen::VectorXcd u_sol_dir_N = cont_space.Interpolate_helmholtz(sol_D, mesh);
                Eigen::VectorXcd u_sol_neu_N = discont_space.Interpolate_helmholtz(sol_N, mesh);
                Eigen::VectorXcd u_sol_N(2*numpanels);
                u_sol_N << u_sol_dir_N, u_sol_neu_N;
                // Computing V matrix
                std::cout << "starting computation of operators" << std::endl;
                Eigen::MatrixXcd M =
                        mass_matrix::GalerkinMatrix(mesh,discont_space,cont_space,order,c_o,0.0);
                std::cout << "mass matrix computed" << std::endl;
                //std::cout << M << std::endl;
                Eigen::MatrixXcd K =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, c_o, 0.0);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, c_o, 0.0);
                std::cout << "hypersingular helmholtz computed" << std::endl;

                // Build rhs for solving
                Eigen::VectorXcd rhs = (0.5*M+K.transpose())*u_sol_neu_N;
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(W);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::cout << "_________________" << std::endl;
                std::cout << sol << std::endl;
                std::cout << "*****************" << std::endl;
                std::cout << (u_sol_dir_N) << std::endl;
                std::cout << "_________________" << std::endl;
                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_data.dat", std::ios_base::app);
                filename << (sol - u_sol_dir_N).lpNorm<2>() << " "
                         << 2 * M_PI / numpanels << std::endl;
                return sol;
            }
            Eigen::VectorXcd solve_3(const ParametrizedMesh &mesh,
                                     std::function<complex<double>(double, double)> u_inc,
                                     std::function<complex<double>(double, double)> u_inc_T,
                                     std::function<complex<double>(double, double)> sol_D,
                                     std::function<complex<double>(double, double)> sol_N,
                                     unsigned order, const double c_o, const double c_i) {
                typedef complex<double> complex_t;
                complex_t ii = complex_t(0, 1);
                unsigned int numpanels = mesh.getNumPanels();
                DiscontinuousSpace<0> discont_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                Eigen::VectorXcd u_sol_dir_N = cont_space.Interpolate_helmholtz(sol_D, mesh);
                Eigen::VectorXcd u_sol_neu_N = discont_space.Interpolate_helmholtz(sol_N, mesh);
                Eigen::VectorXcd u_sol_N(2*numpanels);
                u_sol_N << u_sol_dir_N, u_sol_neu_N;
                // Computing V matrix
                std::cout << "starting computation of operators" << std::endl;
                Eigen::MatrixXcd M =
                        mass_matrix::GalerkinMatrix(mesh,discont_space,cont_space,order,c_o,0.0);
                std::cout << "mass matrix computed" << std::endl;
                Eigen::MatrixXcd K =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, c_o, 0.0);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd V =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, c_o, 0.0);
                std::cout << "single layer helmholtz computed" << std::endl;

                // Build rhs for solving
                Eigen::VectorXcd rhs = (0.5*M-K)*u_sol_dir_N;
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(V);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::ofstream filename;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_data.dat", std::ios_base::app);
                filename << (sol - u_sol_neu_N).lpNorm<2>() << " "
                         << 2 * M_PI / numpanels << std::endl;
                return sol;
            }
            Eigen::VectorXcd solve_2(const ParametrizedMesh &mesh,
                                     std::function<complex<double>(double, double)> u_inc,
                                     std::function<complex<double>(double, double)> u_inc_T,
                                     std::function<complex<double>(double, double)> sol_D,
                                     std::function<complex<double>(double, double)> sol_N,
                                     unsigned order, const double c_o, const double c_i) {
                typedef complex<double> complex_t;
                complex_t ii = complex_t(0, 1);
                int npanels = mesh.getNumPanels();
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
                Eigen::MatrixXcd K_i(npanels,npanels);
                std::string path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd V_i(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/single_layer_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag ) {
                    V_i(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd W_i(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/hypersingular_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    W_i(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_i_adj(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/adjoint_double_layer_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i_adj(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();

                Eigen::MatrixXcd A_i(K_i.rows() + W_i.rows(), K_i.cols() + V_i.cols());
                A_i.block(0, 0, K_i.rows(), K_i.cols()) =  K_i + 0.5* Eigen::MatrixXcd::Identity(K_i.rows(), K_i.cols());
                A_i.block(0, K_i.cols(), V_i.rows(), V_i.cols()) = -V_i;
                A_i.block(K_i.rows(), 0, W_i.rows(), W_i.cols()) = -W_i;
                A_i.block(K_i.rows(), K_i.cols(), K_i.cols(), K_i.rows()) =
                        -K_i_adj + 0.5 * Eigen::MatrixXcd::Identity(K_i.rows(), K_i.cols());
                Eigen::VectorXcd u_inc_dir = cont_space.Interpolate_helmholtz(u_inc, mesh);
                Eigen::VectorXcd u_inc_n = discont_space.Interpolate_helmholtz(u_inc_T, mesh);
                Eigen::VectorXcd u_sol_dir = cont_space.Interpolate_helmholtz(sol_D, mesh);
                Eigen::VectorXcd u_sol_n = discont_space.Interpolate_helmholtz(sol_N, mesh);
                Eigen::VectorXcd u_inc_N(u_inc_dir.size() + u_inc_n.size());
                Eigen::VectorXcd u_sol_N(u_sol_dir.size() + u_sol_n.size());
                u_inc_N << u_inc_dir, u_inc_n;
                u_sol_N << u_sol_dir, u_sol_n;
                // Build rhs for solving
                std::cout <<  u_inc_dir.size() << " " << u_inc_n.size() << std::endl;
                Eigen::VectorXcd rhs = A_i * u_sol_N;
                std::ofstream filename;
                std::cout << "_________________" << std::endl;
                std::cout << rhs << std::endl;
                std::cout << "*****************" << std::endl;
                std::cout << u_sol_N << std::endl;
                std::cout << "_________________" << std::endl;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_data.dat", std::ios_base::app);
                filename << (2*rhs-u_sol_N).segment(0, npanels).lpNorm<2>() << " " << sqrt(2 * M_PI / npanels)
                         << " " << (rhs).segment(0, npanels).lpNorm<2>() / sqrt(2 * M_PI / npanels) << " "
                         <<
                         u_sol_dir.cwiseAbs().mean() << std::endl;

                std::ofstream filename1;
                filename1.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_data.dat", std::ios_base::app);
                filename1 << (2*rhs-u_sol_N).segment(npanels, npanels).lpNorm<2>() << " "
                          << sqrt(2 * M_PI / npanels)
                          << " "
                          << (rhs).segment(npanels, npanels).lpNorm<2>() / sqrt(2 * M_PI / npanels)
                          << " " <<
                          u_sol_n.cwiseAbs().mean() << std::endl;
                return rhs;
            }
            Eigen::VectorXcd solve_1(const ParametrizedMesh &mesh,
                                     std::function<complex<double>(double, double)> u_inc,
                                     std::function<complex<double>(double, double)> u_inc_T,
                                     std::function<complex<double>(double, double)> sol_D,
                                     std::function<complex<double>(double, double)> sol_N,
                                     unsigned order, const double c_o, const double c_i) {
                typedef complex<double> complex_t;
                complex_t ii = complex_t(0, 1);
                int npanels = mesh.getNumPanels();
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
                Eigen::MatrixXcd K_o(npanels,npanels);
                std::string path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_o_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_o(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_i(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/double_layer_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_i_adj(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/adjoint_double_layer_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_i_adj(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd K_o_adj(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/adjoint_double_layer_o_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    K_o_adj(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd V_o(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/single_layer_o_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    V_o(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real,imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd V_i(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/single_layer_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag ) {
                    V_i(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd W_o(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/hypersingular_o_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag ) {
                    W_o(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real, imag);
                    i++;
                    fp_data >> sign >> sign;
                }
                fp_data.close();
                i = 0;
                Eigen::MatrixXcd W_i(npanels,npanels);
                path = "/home/diegorenner/Uni/Thesis/HelmholtzBEM/hypersingular_i_" + std::to_string(npanels) + ".dat";
                fp_data.open(path);
                sign='+';
                while(fp_data >> real >> imag) {
                    W_i(i/npanels,i%npanels) = complex_t((sign=='-')?-real:real, imag);
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
                Eigen::VectorXcd u_inc_dir = cont_space.Interpolate_helmholtz(u_inc, mesh);
                Eigen::VectorXcd u_inc_n = discont_space.Interpolate_helmholtz(u_inc_T, mesh);
                Eigen::VectorXcd u_sol_dir = cont_space.Interpolate_helmholtz(sol_D, mesh);
                Eigen::VectorXcd u_sol_n = discont_space.Interpolate_helmholtz(sol_N, mesh);
                Eigen::VectorXcd u_inc_N(u_inc_dir.size() + u_inc_n.size());
                Eigen::VectorXcd u_sol_N(u_sol_dir.size() + u_sol_n.size());
                u_inc_N << u_inc_dir, u_inc_n;
                u_sol_N << u_sol_dir, u_sol_n;
                // Build rhs for solving
                std::cout <<  u_inc_dir.size() << " " << u_inc_n.size() << std::endl;
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
                filename << (2.0*sol - u_sol_N).segment(0, npanels).lpNorm<2>() << " " << 2 * M_PI / npanels
                         << " " << (sol - u_sol_N).segment(0, npanels).lpNorm<2>() / sqrt(2 * M_PI / npanels) << " "
                         <<
                         u_sol_dir.cwiseAbs().mean() << std::endl;

                std::ofstream filename1;
                filename1.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_data.dat", std::ios_base::app);
                filename1 << ((2.0*sol.real() - 2.0*ii*sol.imag()) - u_sol_N).segment(npanels, npanels).lpNorm<2>() << " "
                          << 2 * M_PI / npanels
                          << " "
                          << (sol - u_sol_N).segment(npanels, npanels).lpNorm<2>() / sqrt(2 * M_PI / npanels)
                          << " " <<
                          u_sol_n.cwiseAbs().mean() << std::endl;
                //std::cout << test.segment(0,u_inc_dir.size()).lpNorm<Eigen::Infinity>() << std::endl;
                return sol;
            }
            Eigen::VectorXcd solve(const ParametrizedMesh &mesh,
                                   std::function<complex<double>(double, double)> u_inc,
                                   std::function<complex<double>(double, double)> u_inc_T,
                                   std::function<complex<double>(double, double)> sol_D,
                                   std::function<complex<double>(double, double)> sol_N,
                                   unsigned order, const double c_o, const double c_i) {
                typedef complex<double> complex_t;
                complex_t ii = complex_t(0, 1);
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd K_o =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, c_o, 0.0);
                std::cout << "double layer helmholtz rhs computed" << std::endl;
                //std::cout << K_o << std::endl;
                Eigen::MatrixXcd K_i =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, c_i, 0.0);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_i =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, c_i, 0.0);
                std::cout << "hypersingular helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_o =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, c_o, 0.0);
                std::cout << "hypersingular helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_o =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, c_o, 0.0);
                std::cout << "single layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_i =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, c_i, 0.0);
                std::cout << "single layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd M =
                        mass_matrix::GalerkinMatrix(mesh,discont_space,cont_space,order,c_o,0.0);
                std::cout << "mass matrix computed" << std::endl;

                Eigen::MatrixXcd A(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A.block(0, 0, K_o.rows(), K_o.cols()) = - K_o + K_i+M;
                A.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o-V_i;
                A.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o-W_i;
                A.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        (K_o-K_i).transpose()+M;
                         ;
                Eigen::MatrixXcd A_o(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
                A_o.block(0, 0, K_o.rows(), K_o.cols()) = -K_o + 0.5*M;
                A_o.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o;
                A_o.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o;
                A_o.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                        K_o.transpose()+0.5*M;

                Eigen::VectorXcd u_inc_dir = cont_space.Interpolate_helmholtz(u_inc, mesh);
                Eigen::VectorXcd u_inc_n = discont_space.Interpolate_helmholtz(u_inc_T, mesh);
                Eigen::VectorXcd u_sol_dir = cont_space.Interpolate_helmholtz(sol_D, mesh);
                Eigen::VectorXcd u_sol_n = discont_space.Interpolate_helmholtz(sol_N, mesh);
                Eigen::VectorXcd u_inc_N(u_inc_dir.size() + u_inc_n.size());
                Eigen::VectorXcd u_sol_N(u_sol_dir.size() + u_sol_n.size());
                u_inc_N << u_inc_dir, u_inc_n;
                u_sol_N << u_sol_dir, u_sol_n;
                // Build rhs for solving
                std::cout << A_o.size() << " " << u_inc_N.size() << std::endl;
                Eigen::VectorXcd rhs = A_o * u_inc_N;
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(A);
                Eigen::VectorXcd sol = dec.solve(rhs);
                std::ofstream filename;
                std::cout << "_________________" << std::endl;
                std::cout << sol << std::endl;
                std::cout << "*****************" << std::endl;
                std::cout << (u_sol_N) << std::endl;
                std::cout << "_________________" << std::endl;
                filename.open("/home/diegorenner/Uni/Thesis/matlab_plots/dirichlet_data_transm.dat", std::ios_base::app);
                filename << (sol - u_sol_N).segment(0, numpanels).lpNorm<2>() << " " << 2 * M_PI / numpanels << std::endl;

                std::ofstream filename1;
                filename1.open("/home/diegorenner/Uni/Thesis/matlab_plots/neumann_data_transm.dat", std::ios_base::app);
                filename1 << (sol.real() -  ii*sol.imag() - u_sol_N).segment(numpanels, numpanels).lpNorm<2>() << " " << 2 * M_PI / numpanels << std::endl;
                return sol;
            }
        } // namespace direct_second_kind

/**
 * This namespace contains the solver using the indirect first kind method which
 * has the variational formulation as given in \f$\eqref{eq:iddirVv}\f$. The
 * Solver uses the lowest order BEM spaces for computation.
 */
    } // namespace parametricbem2d
}
#endif // DIRICHLETHPP
