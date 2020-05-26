#include <complex>
#include <iostream>
#include <fstream>
#include "parametrized_circular_arc.hpp"
#include "solvers.hpp"
#include "generate_solution.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "mass_matrix.hpp"

typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);
int main(int argc, char** argv) {

    // define radius of circle
    double eps = atof(argv[1]);
    std::cout << argv[7] << std::endl;
    int l = atoi(argv[2]);
    double a_n[2*l+1];
    double c_i = atof(argv[3]);
    double c_o = atof(argv[4]);
    double k = atof(argv[5]);//2.75679178324354;
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),eps,0,2*M_PI);
    unsigned order = atoi(argv[6]);
    for( int i = 0; i<2*l+1; i++) {
        a_n[i] = 1./((k*k*(c_o-c_i))*sqrt((2*l+1)*M_PI*eps*eps*(jn(i-l,k)*jn(i-l,k)-jn(i-l-1,k)*jn(i-l+1,k))));
        std::cout << a_n[i] << std::endl;
    }
    std::ofstream filename;
    filename.open(argv[7], std::ofstream::out | std::ofstream::trunc);
    filename.close();
    unsigned n_runs = 7;
    int numpanels[n_runs];
    numpanels[0] = 50;
    for (int i=1; i<n_runs; i++){
        numpanels[i] = 2*numpanels[i-1];
    }
    auto u_i_dir = [&] (double x1, double x2) {
        return sol::u_i(x1, x2, l, eps, a_n, k, c_i);
    };
    auto u_t_dir = [&] (double x1, double x2) {
        return sol::u_t(x1, x2, l, eps, a_n, k, c_i);
    };
    auto u_i_neu = [&] (double x1, double x2) {
        return sol::u_i_N(x1, x2, l, eps, a_n, k, c_i);
    };
    auto u_t_neu = [&] (double x1, double x2) {
        return sol::u_t_N(x1, x2, l, eps, a_n, k, c_i);
    };
    DiscontinuousSpace<0> discont_space;
    ContinuousSpace<1> cont_space;
    // Loop over number of panels
    for (unsigned i = 0; i <= n_runs; i++) {
        ParametrizedMesh mesh(curve.split(numpanels[i]));
        Eigen::VectorXcd sol = tsp::direct_second_kind::solve(
                mesh, u_i_dir, u_i_neu, u_t_dir, u_t_neu, order, k, c_o, c_i);
        Eigen::MatrixXcd M_cont = mass_matrix::GalerkinMatrix(mesh,cont_space,cont_space,order);
        Eigen::MatrixXcd M_discont = mass_matrix::GalerkinMatrix(mesh,discont_space,discont_space,order);
        Eigen::MatrixXcd M(2*numpanels[i],2*numpanels[i]);
        M.block(0,0,numpanels[i],numpanels[i]) = M_cont;
        M.block(numpanels[i],numpanels[i],numpanels[i],numpanels[i]) = M_discont;

        Eigen::VectorXcd u_t_dir_N = cont_space.Interpolate_helmholtz(u_t_dir,mesh);
        Eigen::VectorXcd u_t_neu_N = cont_space.Interpolate_helmholtz(u_t_neu,mesh);
        Eigen::VectorXcd u_t_N(2*numpanels[i]);
        u_t_N << u_t_dir_N, u_t_neu_N;
        filename.open(argv[7], std::ios_base::app);
        filename << mesh.getPanels()[0]->length() << " " << sqrt(abs((sol-u_t_N).dot(M*(sol-u_t_N)))) << std::endl;
        filename.close();

    }
    return 0;
}
