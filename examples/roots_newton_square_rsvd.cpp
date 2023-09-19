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
//#include "singular_values_arnoldi.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"
#include "randsvd.hpp"
#include "find_roots.hpp"
#include "continuous_space.hpp"

#define REFINE
#define PARALLELIZE

// define shorthand for time benchmarking tools, complex data type and immaginary unit
using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

int main(int argc, char** argv) {

    //Eigen::setNbThreads(1);

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

    int nc = 10;
    int nr = 2 * numpanels;
    Eigen::MatrixXcd W = randomized_svd::randGaussian(nr, nc);

    double k_max = 10.0, k_step = (k_max - k_min) / (n_points_k - 1);

    // Inform user of started computation.
#ifdef CMDL
    std::cout << "----------------------------------------------------------------" << std::endl;
    std::cout << "Finding resonances using rSVD approximation and Newton's method." << std::endl;
    std::cout << "Computing on user-defined problem using square domain." << std::endl;
    std::cout << std::endl;
#endif

#ifdef PARALLELIZE
    auto policy = std::execution::par;
    bool profiling = false;
#else
    auto policy = std::execution::seq;
    unsigned total_rsvd_time = 0;
    bool profiling = true;
#endif

    // create objects for assembling solutions operator and its derivatives
    ContinuousSpace<1> cont_space;
    SolutionsOperator so(mesh, order, cont_space, cont_space, profiling);
    GalerkinMatrixBuilder builder(mesh, cont_space, cont_space, order);

    auto tic = high_resolution_clock::now();
#ifdef CMDL
    std::cout << "Bracketing local extrema..." << std::endl;
#endif

    std::vector<size_t> ind(n_points_k);
    std::vector<double> rsv(n_points_k), loc_min, pos_left, val, val_left, val_right;
    std::iota(ind.begin(), ind.end(), 0);

    // Sweep the k interval with subdivision of size n_points_k, do this in parallel.
    // For each value k, approximate the smallest singular value by using rSVD.
    std::transform(policy, ind.cbegin(), ind.cend(), rsv.begin(), [&](size_t i) {
        Eigen::MatrixXcd T;
#ifndef PARALLELIZE
        so.gen_sol_op(builder, k_min + k_step * i, c_o, c_i, T);
        auto tic = high_resolution_clock::now();
#else
        so.gen_sol_op(k_min + k_step * i, c_o, c_i, T);
#endif
        auto res = randomized_svd::sv(T, W, q);
#ifndef PARALLELIZE
        auto toc = high_resolution_clock::now();
        total_rsvd_time += duration_cast<milliseconds>(toc - tic).count();
#endif
        return res;
    });

    // Bracket the local minima of the rSVD curve.
    // If n_points_k is not too large, the obtained intervals will
    // contain the true minima as well (rSVD approximates them
    // with a relative error of about 1e-3).
    // However, if n_points_k is too small, some minima may
    // be missed or the curve may not be convex in the intervals.
    for (size_t i = 0; i < rsv.size() - 2; ++i) {
        double &c = rsv[i+1], L = rsv[i], R = rsv[i+2];
        double k = k_min + i * k_step;
        if (L - c > 0. && R - c > 0.) { // local minimum
            pos_left.push_back(k);
            val.push_back(c);
            val_left.push_back(L);
            val_right.push_back(R);
        }
    }
    size_t disc = 0, loc_min_count = pos_left.size(); // number of local minima
#ifdef CMDL
    std::cout << "Found " << loc_min_count << " candidates. Validating..." << std::endl;
#endif

    auto sv_der = [&](double k) {
        Eigen::MatrixXcd T, T_der;
#ifndef PARALLELIZE
        so.gen_sol_op_1st_der(builder, k, c_o, c_i, T, T_der);
        auto tic = high_resolution_clock::now();
#else
        so.gen_sol_op_1st_der(k, c_o, c_i, T, T_der);
#endif
        auto res = randomized_svd::sv_der(T, T_der, W, q)(1);
#ifndef PARALLELIZE
        auto toc = high_resolution_clock::now();
        total_rsvd_time += duration_cast<milliseconds>(toc - tic).count();
#endif
        return res;
    };
    auto sv_der2 = [&](double k) {
        Eigen::MatrixXcd T, T_der, T_der2;
        Eigen::MatrixXd res(1, 2);
#ifndef PARALLELIZE
        so.gen_sol_op_2nd_der(builder, k, c_o, c_i, T, T_der, T_der2);
        auto tic = high_resolution_clock::now();
#else
        so.gen_sol_op_2nd_der(k, c_o, c_i, T, T_der, T_der2);
#endif
        auto v = randomized_svd::sv_der2(T, T_der, T_der2, W, q);
#ifndef PARALLELIZE
        auto toc = high_resolution_clock::now();
        total_rsvd_time += duration_cast<milliseconds>(toc - tic).count();
#endif
        res(0, 0) = v(1);
        res(0, 1) = v(2);
        return res;
    };

    // Discard candidates which do not approximate local minima
    ind.resize(loc_min_count);
    std::iota(ind.begin(), ind.end(), 0);
    std::vector<bool> accept(loc_min_count, false);
    std::vector<double> der_left(loc_min_count), der_right(loc_min_count);
    auto discard_candidates = [&]() {
        for (size_t i = loc_min_count; i--> 0;) {
            if (!accept[i]) {
                pos_left.erase(pos_left.begin() + i);
                der_left.erase(der_left.begin() + i);
                der_right.erase(der_right.begin() + i);
                val.erase(val.begin() + i);
                val_left.erase(val_left.begin() + i);
                val_right.erase(val_right.begin() + i);
                loc_min_count--;
                disc++;
            }
        }
    };
    std::for_each(policy, ind.cbegin(), ind.cend(), [&](size_t i) {
        der_left[i] = sv_der(pos_left[i]);
        der_right[i] = sv_der(pos_left[i] + 2 * k_step);
        accept[i] = der_left[i] < 0. && der_right[i] > 0.;
    });
    if (disc > 0) {
        discard_candidates();
#ifdef CMDL
        std::cout << "Discarded " << disc << " candidate(s)." << std::endl;
#endif
    }
    loc_min.resize(loc_min_count);
    ind.resize(loc_min_count);

#ifdef REFINE
    // Refine approximations of local minima by using the interpolated quartic
#ifdef CMDL
    std::cout << "Improving candidates..." << std::endl;
#endif
    std::transform(ind.cbegin(), ind.cend(), loc_min.begin(), [&](size_t i) {
        int status = 0;
        double arg, h = k_step, a = pos_left[i], b = pos_left[i] + 2 * h;
        double f1 = val_left[i], f2 = val_right[i], f0 = val[i];
        double d1 = der_left[i], d2 = der_right[i];
        double scale = 0.25 * pow(h, -4), a2 = pow(a, 2);
        double p1 = a2 * (h * (d2 - d1) - 2. * (f1 + f2 - 2. * f0));
        double p2 = a * h * (h * (3. * d2 - 5. * d1) - (9. * f1 + 7. * f2 - 16. * f0));
        double p3 = pow(h, 2) * (2. * h * (d2 - 4. * d1) - (11. * f1 + 5. * f2 - 16. * f0));
        double A = scale * p1 / a2;
        double B = -scale * (4. * p1 + p2) / a;
        double C = scale * (6. * p1 + 3 * p2 + p3);
        double D = d1 - scale * a * (4. * p1 + 3. * p2 + 2. * p3);
        double E = f1 - a * d1 + scale * a2 * (p1 + p2 + p3);
        BrentMinimizer br_min(a, b, acc);
        arg = br_min.local_min_rc(status, 0.);
        while (status) {
            double p = (((A * arg + B) * arg + C) * arg + D) * arg + E;
            arg = br_min.local_min_rc(status, p);
        }
        return arg;
    });
#endif

#ifdef CMDL
    std::cout //<< "Found approximate locations of " << loc_min_count << " local minima." << std::endl
              << "Starting local search..." << std::endl;
#endif
    // Search for minima in the bracketed regions by using Newton-Raphson method
    disc = 0;
    unsigned ict = 0;
    accept.resize(loc_min_count);
    std::for_each(policy, ind.cbegin(), ind.cend(), [&](size_t i) {
        unsigned ic;
        bool rf = false;
        auto fn = [&](double x) {
            if (x == x)
                return double(NAN);
#ifdef REFINE
            return loc_min[i];
#else
            return pos_left[i] + k_step;
#endif
        };
        loc_min[i] = rtsafe(fn, sv_der2, pos_left[i], pos_left[i] + 2 * k_step, acc, rf, ic);
        ict += ic;
        accept[i] = rf;
        if (!rf) ++disc;
    });
    if (disc > 0) {
        discard_candidates();
#ifdef CMDL
        std::cout << "Discarded " << disc << " candidate(s)." << std::endl;
#endif
    }
#ifdef CMDL
    std::cout << "Total Newton iterations taken: " << ict << std::endl;
#endif
    // output results to file
    file_out.open(file_minimas, std::ios_base::app);
    file_out << std::setprecision((int)std::ceil(-std::log10(acc)));
    loc_min_count = loc_min.size();
    for (size_t i = 0; i < loc_min_count; ++i)
        file_out << loc_min[i] << std::endl;
    file_out.close();
    auto toc = high_resolution_clock::now();
#ifdef CMDL
    std::cout << "Total time: " << 1e-3 * duration_cast<milliseconds>(toc - tic).count() << " sec" << std::endl;
#ifndef PARALLELIZE
    std::cout << "Total RSVD time: " << 1e-3 * total_rsvd_time << " sec" << std::endl;
#endif
    std::cout << std::endl;
#endif
    return 0;
}
