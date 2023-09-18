/**
 * \file roots_brent_square_rsvd.cpp
 * \brief This target builds a script that computes minimas in
 * the smallest singular value of the Galerkin BEM
 * approximated solutions operator for the sedond-kind
 * direct BIEs of the Helmholtz transmission problem
 * using the Van Wijngaarden-Dekker-Brent method.
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
 * The resulting file will contain the boundaries of the interval
 * used to compute the root in the first two columns.
 * Then in the next two columns will be the point and the
 * respective function value. The last column will contain the
 * number of iterations used to find the root.
 * The singular values are computed using the Arnoldi algorithm.
 * The user will be updated through the command line about the
 * progress of the algorithm.
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

// define shorthand for time benchmarking tools, complex data type and immaginary unit
using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

// tolerance when verifying root
double epsilon_ver = 1e-3;
// tolerance when finding root
double epsilon_fin = 1e-6;

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
    std::string base_name = "../data/file_roots_brent_square_rsvd_";
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

    std::vector<int> ind(n_points_k), loc_min_iter;
    std::vector<double> rsv(n_points_k),
                        loc_min, loc_min_approx, bracket_left, bracket_right, val, bnd_left, bnd_right
#ifdef FIND_MIN_RSVD
                       ,rsvd_min;
#else
                        ;
#endif

    // Inform user of started computation.
#ifdef CMDL
    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << "Finding resonances using rSVD approximation and Brent's method." << std::endl;
    std::cout << "Computing on user-defined problem using square domain." << std::endl;
    std::cout << std::endl;
#endif

    // create objects for assembling solutions operator and its derivatives
    ContinuousSpace<1> cont_space;
    SolutionsOperator so(mesh, order, cont_space, cont_space, false);
    GalerkinMatrixBuilder builder(mesh, cont_space, cont_space, order);

    std::iota(ind.begin(), ind.end(), 0);
    auto tic = high_resolution_clock::now();
#ifdef CMDL
    std::cout << "Approximating local extrema with randomized SVD..." << std::endl;
#endif
    // Sweep the k interval with subdivision of size n_points_k, do this in parallel.
    // For each value k, approximate the smallest singular value by using rSVD.
    std::transform(std::execution::par_unseq, ind.cbegin(), ind.cend(), rsv.begin(), [&](int i) {
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
            bracket_right.push_back(k + k_step * 2.0);
        }
    }
    if (bracket_right.size() < bracket_left.size())
        bracket_right.push_back(k_max);
    assert (bracket_left.size() == bracket_right.size());

    auto toc = high_resolution_clock::now();

    size_t loc_min_count = bracket_left.size(); // number of local minima
#ifdef CMDL
    std::cout << "Found " << loc_min_count << " candidates for local minima." << std::endl;
    std::cout << "Elapsed time: " << duration_cast<seconds>(toc - tic).count() << " sec" << std::endl;
    std::cout << std::endl;
#endif
    auto sv_eval = [&](double k_in) {
        Eigen::MatrixXcd T_in;
        so.gen_sol_op(builder, k_in, c_o, c_i, T_in);
        return arnoldi::sv(T_in, 1, acc)(0);
    };
    auto sv_eval_der = [&](double k_in) {
        Eigen::MatrixXcd T_in, T_der_in;
        so.gen_sol_op_1st_der(builder, k_in, c_o, c_i, T_in, T_der_in);
        return arnoldi::sv_1st_der(T_in, T_der_in, 1, acc)(0, 1);
    };
    // Search for true minima in the bracketed regions by using
    // Arnoldi iterations with the ZBRENT routine.
    // The builder created at the beginning is used for solutions operator
    // assembly so that it is not created and destroyed in each pass.
    for (size_t i = 0; i < loc_min_count; ++i) {
        double a = bracket_left[i], b = bracket_right[i], r;
        unsigned ic;
        bool rf = false;
#ifdef CMDL
        std::cout << "Local search in [" << a << ", " << b << "]...\t"; std::flush(std::cout);
#endif
        r = zbrent(sv_eval_der, a, b, epsilon_fin, rf, ic);
        if (rf) {
#ifdef CMDL
            std::cout << "Root found at " << r << " in " << ic << " iterations." << std::endl;
#endif
            bnd_left.push_back(a);
            bnd_right.push_back(b);
            loc_min.push_back(r);
            val.push_back(sv_eval(r));
            loc_min_iter.push_back(ic);
        }
    }
#ifdef CMDL
    std::cout << std::endl;
    std::cout << "Found " << loc_min.size() << " local minima." << std::endl;
#endif
    // output minima information to file
    file_out.open(file_minimas, std::ios_base::app);
    loc_min_count = loc_min.size();
    for (size_t i = 0; i < loc_min_count; ++i) {
        file_out << bnd_left[i] << "\t" << bnd_right[i] << "\t"
                 << loc_min[i] << "\t" << val[i] << "\t"
                 << loc_min_iter[i] << std::endl;
    }
    file_out.close();
    toc = high_resolution_clock::now();
#ifdef CMDL
    std::cout << "Total time: " << duration_cast<seconds>(toc - tic).count() << " sec" << std::endl;
    std::cout << std::endl;
#endif
    return 0;
}
