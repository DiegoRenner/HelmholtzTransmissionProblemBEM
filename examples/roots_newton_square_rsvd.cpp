/**
 * \file roots_newton_square_rsvd.cpp
 * \brief This target builds a script that computes minimas in
 * the smallest singular value of the Galerkin BEM
 * approximated solutions operator for the sedond-kind
 * direct BIEs of the Helmholtz transmission problem
 * using the Newton-Raphson method.
 * The scatterer is set to be a square. The results are
 * written to the <tt>data</tt> directory.
 * The script can be run as follows:
 *
 * <tt>
 *  /path/to/roots_brent_square_rsvd \<half side length of square\> \<refraction inside\>
 *     \<refraction outside\> \<initial wavenumber\> \<\#grid points for root search\>
 *     \<\#panels\> \<order of quadrature rule\> \<accuracy of Arnoldi algorithm\>
 *     \<number of subspace iterations\>.
 * </tt>
 *
 * The resulting file will contain the roots in a single column.
 * The user will be updated through the command line about the
 * progress of the algorithm if <tt>CMDL</tt> is set.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 *
 * (c) 2023 Luka Marohnić
 */

#include <complex>
#include <iostream>
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
complex_t ii = complex_t(0,1.);

// tolerance when finding root
double epsilon_fin = 1e-6;
int MAX_CORR = 1 + (int)std::ceil(log2(1.0 / epsilon_fin));

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

    // generate output filename with set parameters
    std::string base_name = "../data/file_roots_newton_square_rsvd_";
    std::string file_plot = base_name + "plot.m";
    std::string suffix = ".dat";
    std::string sep = "_";
    std::string file_minimas = base_name;
    file_minimas.append(argv[2]).append(sep).append(argv[5]).append(sep).append(argv[8]).append(sep).append(argv[9]);
    file_minimas += suffix;
    // clear existing files
    std::ofstream file_out;
    file_out.open(file_minimas, std::ofstream::out | std::ofstream::trunc);
    file_out.close();
    //file_out.open(file_plot, std::ofstream::out | std::ofstream::trunc);
    //file_out.close();

    int nc = 2;
    int nr = 2 * numpanels;
    Eigen::MatrixXcd W = randomized_svd::randGaussian(nr, nc);

    double k_max = 10.0, k_step = (k_max - k_min) / (n_points_k - 1);

    std::vector<size_t> ind(n_points_k), loc_min_pos;
    std::vector<double> rsv(n_points_k), loc_min, bracket_left, bracket_right;

    // Inform user of started computation.
#ifdef CMDL
    std::cout << "----------------------------------------------------------------" << std::endl;
    std::cout << "Finding resonances using rSVD approximation and Newton's method." << std::endl;
    std::cout << "Computing on user-defined problem using square domain." << std::endl;
    std::cout << std::endl;
#endif

    // create objects for assembling solutions operator and its derivatives
    ContinuousSpace<1> cont_space;
    SolutionsOperator so(mesh, order, cont_space, cont_space);
    GalerkinMatrixBuilder builder(mesh, cont_space, cont_space, so.get_GaussQR(), so.get_CGaussQR());

    std::iota(ind.begin(), ind.end(), 0);
    auto tic = high_resolution_clock::now();
#ifdef CMDL
    std::cout << "Approximating local extrema with randomized SVD..." << std::endl;
#endif

    // Sweep the k interval with subdivision of size n_points_k, do this in parallel.
    // For each value k, approximate the smallest singular value by using rSVD.
    std::transform(std::execution::par, ind.cbegin(), ind.cend(), rsv.begin(), [&](size_t i) {
        Eigen::MatrixXcd T;
        so.gen_sol_op(k_min + k_step * i, c_o, c_i, T);
        return randomized_svd::sv(T, W, q);
    });

    // Bracket the local minima of the rSVD curve.
    // If n_points_k is not too large, the obtained intervals will
    // contain the true minima as well (rSVD approximates them
    // with a relative error of about 1e-3).
    // However, if n_points_k is too small, some minima may
    // be missed or the curve may not be convex in the intervals.
    for (size_t i = 0; i < rsv.size() - 2; ++i) {
        double &c = rsv[i+1], L = c - rsv[i], R = rsv[i+2] - c;
        double k = k_min + i * k_step;
        if (L < 0. && R > 0.) { // local minimum
            bracket_left.push_back(k);
            loc_min_pos.push_back(i + 1);
            bracket_right.push_back(k + k_step * 2.0);
        }
    }
    size_t disc = 0, loc_min_count = bracket_left.size(); // number of local minima
#ifdef CMDL
    std::cout << "Found " << loc_min_count << " candidates. Validating..." << std::endl;
#endif

    // Discard candidates which do not approximate local minima
    unsigned N = 2 * numpanels;
    Eigen::MatrixXcd Tall(N, N * loc_min_count), Tall_right(N, N * loc_min_count);
    Eigen::MatrixXcd Tall_der(N, N * loc_min_count), Tall_der_right(N, N * loc_min_count);
    ind.resize(loc_min_count);
    std::iota(ind.begin(), ind.end(), 0);
    std::for_each(std::execution::par, ind.cbegin(), ind.cend(), [&](size_t i) {
        Eigen::MatrixXcd T_left, T_right, T_der_left, T_der_right;
        so.gen_sol_op_1st_der(bracket_left[i], c_o, c_i, T_left, T_der_left);
        so.gen_sol_op_1st_der(bracket_right[i], c_o, c_i, T_right, T_der_right);
        Tall.block(0, i * N, N, N) = T_left;
        Tall_right.block(0, i * N, N, N) = T_right;
        Tall_der.block(0, i * N, N, N) = T_der_left;
        Tall_der_right.block(0, i * N, N, N) = T_der_right;
    });
    for (size_t i = loc_min_count; i--> 0;) {
        const Eigen::MatrixXcd &T = Tall.block(0, i * N, N, N), &T_der = Tall_der.block(0, i * N, N, N);
        const Eigen::MatrixXcd &T_right = Tall_right.block(0, i * N, N, N), &T_der_right = Tall_der_right.block(0, i * N, N, N);
        if (arnoldi::sv_1st_der(T, T_der, 1, acc)(0, 1) * arnoldi::sv_1st_der(T_right, T_der_right, 1, acc)(0, 1) > 0.) {
            bracket_left.erase(bracket_left.begin() + i);
            bracket_right.erase(bracket_right.begin() + i);
            loc_min_pos.erase(loc_min_pos.begin() + i);
            loc_min_count--;
            disc++;
        }
    }
#ifdef CMDL
    if (disc > 0)
        std::cout << "Discarded " << disc << " candidates." << std::endl;
#endif

    loc_min.resize(loc_min_count);
#ifdef BRENT_REFINE
    // Find approximations of local minima by using Brent's algorithm.
#ifdef CMDL
    std::cout << "Refining with Brent's method..." << std::endl;
#endif
    ind.resize(loc_min_count);
    auto rsvd_sv = [&](double k) {
        Eigen::MatrixXcd T;
        so.gen_sol_op(builder, k, c_o, c_i, T);
        return randomized_svd::sv(T, W, q);
    };
    std::transform(ind.cbegin(), ind.cend(), loc_min.begin(), [&](size_t i) {
        int status = 0;
        double arg, a = bracket_left[i], b = bracket_right[i];
        BrentMinimizer brent_minimizer(a, b, epsilon_fin);
        arg = brent_minimizer.local_min_rc(status, 0.);
        while (status)
            arg = brent_minimizer.local_min_rc(status, rsvd_sv(arg));
        return arg;
    });
#else
    for (size_t i = 0; i < loc_min_count; ++i) {
        size_t p = loc_min_pos[i];
        loc_min[i] = k_min + p * k_step;
#if 0
        if (p > 1 && p < n_points_k - 1) {
            double num = -rsv[p - 2] + 8. * rsv[p - 1] - 8. * rsv[p + 1] + rsv[p + 2];
            double den = -rsv[p - 2] + 16. * rsv[p - 1] -30. * rsv[p] + 16. * rsv[p + 1] - rsv[p + 2];
            double pos = loc_min[i] - k_step * num / den;
            if (std::abs(loc_min[i] - pos) < k_step)
                loc_min[i] = pos;
        }
#endif
    }
#endif

#ifdef CMDL
    std::cout << "Found approximate locations of " << loc_min_count << " local minima." << std::endl;
#endif
    // Search for true minima in the bracketed regions by using
    // Arnoldi iterations with Newton-Raphson method.
    // The builder created at the beginning is used for solutions operator
    // assembly so that it is not created and destroyed in each pass.
    int ic = 1 + (int)std::round(-log10(epsilon_fin)/2. - 1.); // number of Newton iterations
#ifdef CMDL
    std::cout << "Refining with Newton's method (" << ic << " iterations)..." << std::endl;
#endif
    Eigen::MatrixXcd Tall_der2 (N, N * loc_min_count);
    disc = 0;
    for (size_t j = 0; j < ic; ++j) {
        std::for_each(std::execution::par, ind.cbegin(), ind.cend(), [&](size_t i) {
            Eigen::MatrixXcd T, T_der, T_der2;
            so.gen_sol_op_2nd_der(loc_min[i], c_o, c_i, T, T_der, T_der2);
            Tall.block(0, i * N, N, N) = T;
            Tall_der.block(0, i * N, N, N) = T_der;
            Tall_der2.block(0, i * N, N, N) = T_der2;
        });
        for (size_t i = loc_min_count; i--> 0;) {
            Eigen::MatrixXcd T = Tall.block(0, i * N, N, N);
            Eigen::MatrixXcd T_der = Tall_der.block(0, i * N, N, N);
            Eigen::MatrixXcd T_der2 = Tall_der2.block(0, i * N, N, N);
            double a = bracket_left[i], b = bracket_right[i], update = NAN, d1, d2;
            if (b - a < epsilon_fin)
                continue; // refinement done by bisection
            size_t corr_count = 0;
            while (true) {
                auto der = arnoldi::sv_2nd_der(T, T_der, T_der2, 1, acc).block(0, 1, 1, 2);
                if (update == update) {
                    if (der(0, 0) < 0. && d1 > 0.)
                        a = update, b = loc_min[i];
                    else if (der(0, 0) > 0. && d1 < 0.)
                        a = loc_min[i], b = update;
                    else if (d1 > 0.)
                        b = update;
                    else a = update;
                    bracket_left[i] = a;
                    bracket_right[i] = b;
                    if (b - a < epsilon_fin) {
                        loc_min[i] = (a + b) / 2.;
                        break;
                    }
                    loc_min[i] = update;
                }
                d1 = der(0, 0);
                d2 = der(0, 1);
                double k = loc_min[i] - d1 / d2;
                if (k > b || k < a) {
                    ++corr_count;
                    if (corr_count > MAX_CORR)
                        break; // discard this candidate
                    // bisect and recompute solutions matrix and derivatives,
                    // then do Newton step again
                    if (d1 < 0.)
                        update = (loc_min[i] + b) / 2.;
                    else update = (loc_min[i] + a) / 2.;
                    so.gen_sol_op_2nd_der(builder, update, c_o, c_i, T, T_der, T_der2);
                } else {
                    // Newton step accepted
                    loc_min[i] = k;
                    break;
                }
            }
            if (corr_count > MAX_CORR) {
                // discard candidate
                loc_min.erase(loc_min.begin() + i);
                ind.erase(ind.begin() + i);
                bracket_left.erase(bracket_left.begin() + i);
                bracket_right.erase(bracket_right.begin() + i);
                loc_min_count--;
                disc++;
            }
        }
    }
#ifdef CMDL
    if (disc > 0)
        std::cout << "Discarded " << disc << " candidates." << std::endl;
#endif
    // output results to file
    file_out.open(file_minimas, std::ios_base::app);
    loc_min_count = loc_min.size();
    for (size_t i = 0; i < loc_min_count; ++i)
        file_out << loc_min[i] << std::endl;
    file_out.close();
    auto toc = high_resolution_clock::now();
#ifdef CMDL
    std::cout << "Total time: " << duration_cast<seconds>(toc - tic).count() << " sec" << std::endl;
    std::cout << std::endl;
#endif
    return 0;
}
