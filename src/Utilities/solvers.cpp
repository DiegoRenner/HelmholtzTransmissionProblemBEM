/**
 * \file
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem of the form given in \f$\eqref{eq:dirbvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */
#include "solvers.hpp"
#include "mass_matrix.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "double_layer.hpp"
#include "hypersingular.hpp"
#include "parametrized_mesh.hpp"
#include "single_layer.hpp"
#include <fstream>
#include <iostream>

    typedef std::complex<double> complex_t;
    namespace bvp {
        namespace direct_first_kind {
            Eigen::VectorXcd solve_neumann(const ParametrizedMesh &mesh,
                                             const std::function<complex_t(double, double)> u_dir,
                                             const std::function<complex_t(double, double)> u_neu,
                                             const unsigned int order,
                                             const double k,
                                             const double c) {
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
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k, c);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c);
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
                                           const std::function<complex_t(double, double)> u_dir,
                                           const std::function<complex_t(double, double)> u_neu,
                                           const unsigned order,
                                           const double k,
                                           const double c) {
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
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k, c);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd V =
                        single_layer_helmholtz::GalerkinMatrix(mesh,discont_space,order,k, c);
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
            Eigen::VectorXcd solve(const ParametrizedMesh &mesh,
                                   const std::function<complex_t(double, double)> u_inc_dir,
                                   const std::function<complex_t(double, double)> u_inc_neu,
                                   const std::function<complex_t(double, double)> sol_dir,
                                   const std::function<complex_t(double, double)> sol_neu,
                                   const unsigned order,
                                   const double k,
                                   const double c_o,
                                   const double c_i) {
                int numpanels = mesh.getNumPanels();
                // Same trial and test spaces
                DiscontinuousSpace<0> discont_space;
                //DiscontinuousSpace<0> test_space;
                // Space used for interpolation of Dirichlet data
                ContinuousSpace<1> cont_space;
                // Computing V matrix
                std::cout << "starting computation of operators for " << numpanels << " panels." << std::endl;
                Eigen::MatrixXcd K_o =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k, c_o);
                std::cout << "double layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd K_i =
                        double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k, c_i);
                std::cout << "double layer helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_i =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
                std::cout << "hypersingular helmholtz computed" << std::endl;
                Eigen::MatrixXcd W_o =
                        hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k, c_o);
                std::cout << "hypersingular helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_o =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k, c_o);
                std::cout << "single layer helmholtz rhs computed" << std::endl;
                Eigen::MatrixXcd V_i =
                        single_layer_helmholtz::GalerkinMatrix(mesh, discont_space, order, k, c_i);
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
                Eigen::VectorXcd u_inc_N(2*numpanels);
                u_inc_N << u_inc_dir_N, u_inc_neu_N;
                // Build rhs for solving
                Eigen::VectorXcd rhs = (A_o * u_inc_N);
                // Solving for coefficients
                Eigen::HouseholderQR<Eigen::MatrixXcd> dec(A);
                Eigen::VectorXcd sol = dec.solve(rhs);
                //std::cout << "-----------------"<< std::endl;
                //std::cout << sol.segment(0,2*numpanels).transpose() << std::endl;
                //std::cout << "**************************"<< std::endl;
                //std::cout << u_sol_N.segment(0,2*numpanels).transpose() << std::endl;
                //std::cout << "-----------------"<< std::endl;

                //std::cout << "****************************" << std::endl;
                return sol;
            }

        } // namespace direct_second_kind
    } // namespace tsp
