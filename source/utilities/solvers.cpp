#include "solvers.hpp"
#include "mass_matrix.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "parametrized_mesh.hpp"
#include "single_layer.hpp"
#include "double_layer.hpp"
#include "hypersingular.hpp"

typedef std::complex<double> complex_t;
namespace bvp {
    namespace direct_first_kind {
        Eigen::VectorXcd solve_dirichlet(const ParametrizedMesh &mesh,
                                         const std::function<complex_t(double, double)> u_dir,
                                         const unsigned order,
                                         const double k) {
            // define #panels for convenience
            unsigned int numpanels = mesh.getNumPanels();
            // compute FEM-spaces of lowest order
            ContinuousSpace<1> cont_space;
            DiscontinuousSpace<0> discont_space;
            // compute interpolation coefficients for dirichlet data
            Eigen::VectorXcd u_dir_N = cont_space.Interpolate_helmholtz(u_dir, mesh);
            // compute operators for first kind direct Dirichlet problem BIE
            Eigen::MatrixXcd M =
                    mass_matrix::GalerkinMatrix(mesh,discont_space,cont_space,order);
            Eigen::MatrixXcd K =
                    double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k,1);
            Eigen::MatrixXcd V =
                    single_layer_helmholtz::GalerkinMatrix(mesh,discont_space,order,k,1);
            // build rhs for solving
            Eigen::VectorXcd rhs = (0.5*M-K)*u_dir_N;
            // solve for coefficients
            Eigen::HouseholderQR<Eigen::MatrixXcd> dec(-V);
            Eigen::VectorXcd sol = dec.solve(rhs);
            return sol;
        }
        Eigen::VectorXcd solve_neumann(const ParametrizedMesh &mesh,
                                       const std::function<complex_t(double, double)> u_neu,
                                       const unsigned int order,
                                       const double k) {
            // define #panels for convenience
            unsigned int numpanels = mesh.getNumPanels();
            // compute FEM-spaces of lowest order
            ContinuousSpace<1> cont_space;
            DiscontinuousSpace<0> discont_space;
            // compute interpolation coefficients of Neumann data
            Eigen::VectorXcd u_neu_N = discont_space.Interpolate_helmholtz(u_neu, mesh);
            // compute operators for first kind direct Neumann problem BIE
            Eigen::MatrixXcd M =
                    mass_matrix::GalerkinMatrix(mesh,discont_space,cont_space,order);
            Eigen::MatrixXcd K =
                    double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, discont_space, order, k, 1);
            Eigen::MatrixXcd W =
                    hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, 1);
            // build rhs for solving
            Eigen::VectorXcd rhs = (0.5*M+K.transpose())*u_neu_N;
            // solve for coefficients
            Eigen::HouseholderQR<Eigen::MatrixXcd> dec(-W);
            Eigen::VectorXcd sol = dec.solve(rhs);
            return sol;
        }
    } // namespace direct_first_kind
} // namespace bvp
namespace tp {
    namespace direct_second_kind {
        Eigen::VectorXcd solve(const ParametrizedMesh &mesh,
                               const std::function<complex_t(double, double)> u_inc_dir,
                               const std::function<complex_t(double, double)> u_inc_neu,
                               const unsigned order,
                               const double k,
                               const double c_o,
                               const double c_i) {
            // define number of panels for convenience
            int numpanels = mesh.getNumPanels();
            // space used for interpolation of Dirichlet data
            ContinuousSpace<1> cont_space;
            // compute operators of second kind direct BIEs for the Helmholtz Transmission problem
            Eigen::MatrixXcd K_o =
                    double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_o);
            Eigen::MatrixXcd K_i =
                    double_layer_helmholtz::GalerkinMatrix(mesh, cont_space, cont_space, order, k, c_i);
            Eigen::MatrixXcd W_i =
                    hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
            Eigen::MatrixXcd W_o =
                    hypersingular_helmholtz::GalerkinMatrix(mesh, cont_space, order,k, c_o);
            Eigen::MatrixXcd V_o =
                    single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_o);
            Eigen::MatrixXcd V_i =
                    single_layer_helmholtz::GalerkinMatrix(mesh, cont_space, order, k, c_i);
            Eigen::MatrixXcd M_cont =
                    mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
            // Build matrices for solving linear system of equations
            Eigen::MatrixXcd A(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
            A.block(0, 0, K_o.rows(), K_o.cols()) = (-K_o + K_i)+M_cont;
            A.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = (V_o-V_i);
            A.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o-W_i;
            A.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                    (K_o-K_i).transpose()+M_cont;
            // Build matrices for right hand side
            Eigen::MatrixXcd A_o(K_o.rows() + W_o.rows(), K_o.cols() + V_o.cols());
            A_o.block(0, 0, K_o.rows(), K_o.cols()) = -K_o + 0.5*M_cont;
            A_o.block(0, K_o.cols(), V_o.rows(), V_o.cols()) = V_o;
            A_o.block(K_o.rows(), 0, W_o.rows(), W_o.cols()) = W_o;
            A_o.block(K_o.rows(), K_o.cols(), K_o.cols(), K_o.rows()) =
                    K_o.transpose()+0.5*M_cont;

            // Build vectors from incoming wave data for right hand side
            Eigen::VectorXcd u_inc_dir_N = cont_space.Interpolate_helmholtz(u_inc_dir, mesh);
            Eigen::VectorXcd u_inc_neu_N = cont_space.Interpolate_helmholtz(u_inc_neu, mesh);
            Eigen::VectorXcd u_inc_N(2*numpanels);
            u_inc_N << u_inc_dir_N, u_inc_neu_N;
            // compute right hand side
            Eigen::VectorXcd rhs = (A_o * u_inc_N);
            // Solving for coefficients
            Eigen::HouseholderQR<Eigen::MatrixXcd> dec(A);
            Eigen::VectorXcd sol = dec.solve(rhs);
            return sol;
        }
    } // namespace direct_second_kind
} // namespace tsp
