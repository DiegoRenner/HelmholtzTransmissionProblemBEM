/**
 * \file rsvd_test.cpp
 * \brief This target builds a script that compares
 * rSVD and Arnoldi curves.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 *
 * (c) 2023 Luka Marohnić
 */

#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <execution>
#include <algorithm>
#include <string>
#include "parametrized_line.hpp"
#include "singular_values_arnoldi.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "randsvd.hpp"
#include "find_roots.hpp"
#include "continuous_space.hpp"

//#define BRENT_REFINE

// define shorthand for time benchmarking tools, complex data type and immaginary unit
using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0, 1.);

// tolerance when finding root
double epsilon_fin = 1e-12;

int main(int argc, char** argv) {

    // check whether we have the correct number of input arguments
    if (argc < 10)
        throw std::runtime_error("Too few input arguments!");
    if (argc > 10)
        throw std::runtime_error("Too many input arguments!");

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[2]);
    double c_o = atof(argv[3]);
    double k_min = atof(argv[4]);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_k = atoi(argv[5]);
    unsigned numpanels;
    numpanels = atoi(argv[6]);
    // compute mesh for numpanels
    using PanelVector = PanelVector;
    // corner points for the square
    Eigen::RowVectorXd x1(2);
    x1 << 0,0; // point (0,0)
    Eigen::RowVectorXd x2(2);
    x2 << eps, 0; // point (1,0)
    Eigen::RowVectorXd x3(2);
    x3 << eps, eps; // point (1,0.5)
    Eigen::RowVectorXd x4(2);
    x4 << 0, eps; // point (0,1.5)
    // parametrized line segments forming the edges of the polygon
    ParametrizedLine line1(x1, x2);
    ParametrizedLine line2(x2, x3);
    ParametrizedLine line3(x3, x4);
    ParametrizedLine line4(x4, x1);
    // splitting the parametrized lines into panels for a mesh to be used for
    // BEM (Discretization).
    PanelVector line1panels = line1.split(numpanels/4);
    PanelVector line2panels = line2.split(numpanels/4);
    PanelVector line3panels = line3.split(numpanels/4);
    PanelVector line4panels = line4.split(numpanels/4);
    PanelVector panels;
    // storing all the panels in order so that they form a polygon
    panels.insert(panels.end(), line1panels.begin(), line1panels.end());
    panels.insert(panels.end(), line2panels.begin(), line2panels.end());
    panels.insert(panels.end(), line3panels.begin(), line3panels.end());
    panels.insert(panels.end(), line4panels.begin(), line4panels.end());
    // construction of a ParametrizedMesh object from the vector of panels
    ParametrizedMesh mesh(panels);

    // define order of quadrature rule used to compute matrix entries and which singular value to evaluate
    unsigned order = atoi(argv[7]);

    // define accurracy of arnoldi algorithm
    double acc = std::max(atof(argv[8]), std::numeric_limits<double>::epsilon());

    // define the number of subspace iterations
    int q = atoi(argv[9]);

    std::string mfile = "../data/rsvd_test.m";
    std::ofstream file_out;
    file_out.open(mfile, std::ofstream::out | std::ofstream::trunc);
    file_out.close();

    int nc = 4;
    int nr = 2 * numpanels;
    Eigen::MatrixXcd W = randomized_svd::randGaussian(nr, nc);

    double k_max = k_min + .5, k_step = (k_max - k_min) / (n_points_k - 1);

    std::vector<size_t> ind(n_points_k);
    std::vector<double> rsv(n_points_k), asv(n_points_k), rsv_min, asv_min, bracket_left, bracket_right;
    std::vector<double>::const_iterator it;

    // create objects for assembling solutions operator and its derivatives
    ContinuousSpace<1> cont_space;
    SolutionsOperator so(mesh, order, cont_space, cont_space);
    GalerkinMatrixBuilder builder(mesh, cont_space, cont_space, so.get_GaussQR(), so.get_CGaussQR());

    std::iota(ind.begin(), ind.end(), 0);

    Eigen::MatrixXcd Tall(nr, nr * n_points_k);
    auto rsvd_sv = [&](double k) {
        Eigen::MatrixXcd T;
        so.gen_sol_op(k, c_o, c_i, T);
        return randomized_svd::sv(T, W, q);
    };
    auto rsvd_sv_der = [&](double k) {
        Eigen::MatrixXcd T, T_der;
        so.gen_sol_op_1st_der(k, c_o, c_i, T, T_der);
        return randomized_svd::sv_der(T, T_der, W, q)(1);
    };
    auto rsvd_sv_both = [&](double k) {
        Eigen::MatrixXcd T, T_der, T_der2;
        so.gen_sol_op_2nd_der(k, c_o, c_i, T, T_der, T_der2);
        Eigen::MatrixXd res(1, 2);
        auto v = randomized_svd::sv_der2(T, T_der, T_der2, W, q);
        res(0, 0) = v(1);
        res(0, 1) = v(2);
        return res;
    };

#ifdef CMDL
    std::cout << "Computing rSVD and Arnoldi curves for c_i = " << c_i << "..." << std::endl;
#endif
    auto tic = high_resolution_clock::now();
    std::for_each(std::execution::par, ind.cbegin(), ind.cend(), [&](size_t i) {
        Eigen::MatrixXcd T_in;
        so.gen_sol_op(k_min + k_step * i, c_o, c_i, T_in);
        Tall.block(0, i * nr, nr, nr) = T_in;
    });
    auto toc = high_resolution_clock::now();
    auto dur_assembly = duration_cast<milliseconds>(toc - tic);
    auto dur_rsv = duration_cast<milliseconds>(toc - tic);
    std::cout << "Assembly: " << dur_assembly.count() << std::endl;
#if 0
    std::cout << "Hankel time: " << builder.getHankelTime() << ", ratio: "
              << (100 * builder.getHankelTime()) / dur_assembly.count() << "%" << std::endl;
    std::cout << "Interaction matrix assembly time: " << (100 * builder.getInteractionMatrixTime()) / dur_assembly.count() << "%" << std::endl;
    std::cout << "Mapping time: " << (100 * builder.getMappingTime()) / dur_assembly.count() << "%" << std::endl;
#endif
    tic = high_resolution_clock::now();
    for (size_t i = 0; i < n_points_k; ++i) {
        const Eigen::MatrixXcd &T = Tall.block(0, i * nr, nr, nr);
        asv[i] = arnoldi::sv(T, 1, acc)(0);
    }
    toc = high_resolution_clock::now();
    auto dur_asv = duration_cast<milliseconds>(toc - tic);
    tic = high_resolution_clock::now();
    std::for_each(std::execution::par, ind.cbegin(), ind.cend(), [&](size_t i) {
        const Eigen::MatrixXcd &T = Tall.block(0, i * nr, nr, nr);
        rsv[i] = randomized_svd::sv(T, W, q);
    });
    toc = high_resolution_clock::now();
    std::cout << "Arnoldi: " << dur_asv.count() << std::endl;
    std::cout << "RSVD: " << dur_rsv.count() << std::endl;
    std::cout << "Ratio: " << double(dur_assembly.count()) / double(dur_rsv.count()) << std::endl;
    std::vector<double> err(n_points_k);
    for (size_t i = 0; i < n_points_k; ++i) {
        err[i] = std::abs(asv[i] - rsv[i]) / std::abs(asv[i]);
    }
    auto err_sorted = err;
    std::sort(err_sorted.begin(), err_sorted.end());
    std::cout << "Error: min " << err_sorted[0] << ", max " << err_sorted[n_points_k - 1]
              << ", median " << err_sorted[n_points_k / 2] << std::endl;

    file_out.open(mfile, std::ios_base::app);
    file_out << std::setprecision(18);
    file_out << "k = [";
    for (size_t i = 0; i < n_points_k; ++i) {
        if (i > 0)
            file_out << ",";
        file_out << k_min + i * k_step;
    }
    file_out << "];" << std::endl << "rsv = [";
    for (it = rsv.begin(); it != rsv.end(); ++it) {
        if (it != rsv.begin())
            file_out << ",";
        file_out << *it;
    }
    file_out << "];" << std::endl << "asv = [";
    for (it = asv.begin(); it != asv.end(); ++it) {
        if (it != asv.begin())
            file_out << ",";
        file_out << *it;
    }
    file_out << "];" << std::endl << "err = [";
    for (it = err.begin(); it != err.end(); ++it) {
        if (it != err.begin())
            file_out << ",";
        file_out << *it;
    }
    file_out << "];" << std::endl;
    file_out.close();

    return 0;
}
