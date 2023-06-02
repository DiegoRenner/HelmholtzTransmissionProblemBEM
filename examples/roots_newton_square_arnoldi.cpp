/**
 * \file roots_newton_square_arnoldi.cpp
 * \brief This target builds a script that computes minimas in the smallest singular value of the
 * Galerkin BEM approximated solutions operator for the second-kind direct BIEs of the Helmholtz
 * transmission problem using the Newton-Raphson method.
 * The scatterer is set to be a square.
 * The results are written to disk.
 * The script can be run as follows:
 *
 * <tt>
 *  /path/to/roots_newton_circle \<half side length of square\> \<refraction inside\>
 *     \<refraction outside\> \<initial wavenumber\> \<\#grid points for root search\>
 *     \<\#panels\> \<order of quadrature rule\> \<outputfile\>.
 * </tt>
 *
 * The resulting file will contain the left boundary of the 
 * interval used to compute the root in the first column. 
 * Then in the next three columns will be the point, 
 * the function value and the derivative at which the root was found.
 * The last column will contain the number of iterations used to find the root.
 * If no root was found the last four columns will be set to \f$\verb|NAN|\f$.
 * The singular values and their derivatives are computed using the Arnoldi algorithm.
 * The user will be updated through the command line about the
 * progress of the algorithm
 * if \f$ \verb|-DCMDL| \f$ is set.
 *
 * This File is a part of the HelmholtzTransmissionProblemBEM library.
 */

#include <complex>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <execution>
#include <algorithm>
#include <string>
#include <limits>
#include "parametrized_line.hpp"
#include "singular_values_arnoldi.hpp"
#include "find_roots.hpp"
#include "gen_sol_op.hpp"

// define shorthand for time benchmarking tools, complex data type and immaginary unit
using namespace std::chrono;
typedef std::complex<double> complex_t;
complex_t ii = complex_t(0,1.);

// tolerance when verifying root
double epsilon_ver = 1e-3;
// tolerance when finding root
double epsilon_fin = 1e-6;

// Create a standard normal random matrix with NR rows and NC columns
Eigen::MatrixXcd standard_normal_random_matrix(int nr, int nc) {
    std::random_device rd {};
    std::mt19937 gen { rd() };
    std::normal_distribution<> d { 0, 1 };
    Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(nr, nc);
    for (int i = 0; i < nr; ++i) for (int j = 0; j < nc; ++j) {
        W(i, j) = complex_t(d(gen), d(gen));
    }
    return W;
}

// Approximate the smallest singular value by randomized SVD
inline double rsv(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &W) {
    int nr = W.rows(), nc = W.cols();
    Eigen::MatrixXcd U, V, Q, B, C, thinQ = Eigen::MatrixXcd::Identity(nr, nc);
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu_decomp(T * T.transpose().conjugate() * T);
    V = lu_decomp.solve(W);
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr(V);
    Q = qr.householderQ() * thinQ;
    B = T.transpose().conjugate() * (T * lu_decomp.solve(Q));
    C = Q.transpose().conjugate() * B;
    Eigen::BDCSVD<Eigen::MatrixXcd> svd(C * Q.transpose().conjugate());
    return 1.0 / svd.singularValues()(0);
}

int main(int argc, char** argv){

    // define radius of circle refraction index and initial wavenumber
    double eps = atof(argv[1]);
    double c_i = atof(argv[2]);
    double c_o = atof(argv[3]);
    double k_min = atof(argv[4]);

    // define mesh in space and on wavenumber on which to perform verification
    unsigned n_points_k = atoi(argv[5]);
    unsigned numpanels = atoi(argv[6]);
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
    double acc = atof(argv[8]);

    // generate output filename with set parameters
    std::string base_name = "../data/file_roots_newton_square_arnoldi_";
    std::string suffix = ".dat";
    std::string divider = "_";
    std::string file_minimas = base_name.append(argv[2]).append(divider).append(argv[5])
                                       .append(divider).append(argv[8]) + suffix;
    // clear existing file
    std::ofstream file_out;
    file_out.open(file_minimas, std::ofstream::out | std::ofstream::trunc);
    file_out.close();

    int nc = 2;
    int nr = 2 * numpanels;
    Eigen::MatrixXcd W = standard_normal_random_matrix(nr, nc);

    double k_max = 4.0, k_step = (k_max - k_min) / n_points_k;

    std::cout << "Bracketing local minima... [--------------------]";

    auto tic = chrono::high_resolution_clock::now();

    unsigned ksize = 1 + (int)round((k_max - k_min) / k_step), finished = 0;
    std::vector<double> kv(ksize), res(ksize), loc_min, loc_min_approx, bnd(1, k_min);
    std::generate(kv.begin(), kv.end(), [&]() { static double k = k_min - k_step; return k += k_step; });

    std::transform(std::execution::par_unseq, kv.cbegin(), kv.cend(), res.begin(), [&](double k) {
        Eigen::MatrixXcd T = gen_sol_op(mesh, order, k, c_o, c_i);
        double sv = rsv(T, W);
        ++finished;
        size_t progress = (100 * finished) / ksize, l1 = progress / 5, l2 = 20 - l1;
        std::cout << "\rBracketing local minima... [" << std::string(l1, '=') << std::string(l2, '-') << "] (" << progress << "%)";
        std::flush(std::cout);
        return sv;
    });

    std::cout << "\rBracketing local minima... Done.                          " << std::endl;

    for (size_t i = 0; i < res.size() - 2; ++i) {
        double &c = res[i+1], L = c - res[i], R = res[i+2] - c;
        double k = k_min + i * k_step;
        if (L < 0. && R > 0.) { // local minimum
            loc_min_approx.push_back(k);
        }
        if (L > 0. && R < 0.) { // local maximum
            bnd.push_back(k);
        }
    }
    bnd.push_back(k_max);
    auto toc = chrono::high_resolution_clock::now();

    std::cout << "Found approximate locations of " << loc_min_approx.size() << " local minima." << std::endl;
    std::cout << "Elapsed time: " << duration_cast<seconds>(toc - tic).count() << " sec" << std::endl;

    // Routine for computing the smallest singular value by Arnoldi iterations
    auto asv = [&](double k) {
        Eigen::MatrixXcd T_in;
        Eigen::MatrixXcd T_der_in;
        Eigen::MatrixXcd T_der2_in;
        T_in = gen_sol_op(mesh, order, k , c_o, c_i);
        T_der_in = gen_sol_op_1st_der(mesh, order, k , c_o, c_i);
        T_der2_in = gen_sol_op_2nd_der(mesh, order, k , c_o, c_i);
        return arnoldi::sv_2nd_der(T_in, T_der_in, T_der2_in, 1, acc).block(0, 1, 1, 2);
    };

    std::vector<double> d1(loc_min_approx.size()), d2(loc_min_approx.size());
    std::vector<double>::iterator d1t = d1.begin(), d2t = d2.begin();
    std::vector<int> elim_pts;

    int n_changed = 0;

    // refine starting points
    for (std::vector<double>::iterator it = loc_min_approx.begin(); it != loc_min_approx.end(); ++it, ++d1t, ++d2t) {
        bool is_changed = false;
        while (true) {
            Eigen::MatrixXd der = asv(*it);
            *d1t = der(0, 0);
            *d2t = der(0, 1);
            if (*d2t <= 0.) {
                *it += k_step * (*d1t < 0. ? 1.0 : -1.0);
                if (*it <= k_min || *it >= k_max) {
                    elim_pts.push_back(it - loc_min_approx.begin());
                    is_changed = false;
                    break;
                }
                is_changed = true;
            } else break;
        }
        if (is_changed)
            ++n_changed;
    }

    if (n_changed > 0)
        std::cout << "Improved " << n_changed << " starting point(s)" << std::endl;

    if (!elim_pts.empty()) {
        std::cout << "Discarded " << elim_pts.size() << " starting point(s)" << std::endl;
        while (!elim_pts.empty()) {
            loc_min_approx.erase(loc_min_approx.begin() + elim_pts.back());
            d1.erase(loc_min_approx.begin() + elim_pts.back());
            d2.erase(loc_min_approx.begin() + elim_pts.back());
            elim_pts.pop_back();
        }
    }

    d1t = d1.begin();
    d2t = d2.begin();

    for (std::vector<double>::const_iterator it = loc_min_approx.begin(); it != loc_min_approx.end(); ++it, ++d1t, ++d2t) {
        double k0 = *it, k1, lb, ub;
        std::vector<double>::const_iterator bt = std::lower_bound(bnd.begin(), bnd.end(), k0);
        ub = *bt;
        lb = *(--bt);
        size_t iter = 0;
        bool failed = false;
        std::cout << "Finding local minimum " << 1 + int(it - loc_min_approx.begin())
                  << ", k0 = " << *it << ", " << lb << " < k < " << ub << ": ";
        std::flush(std::cout);
        while (true) {
            k1 = k0 - *d1t / *d2t;
            ++iter;
            if (std::abs(k1 - k0) < epsilon_fin)
                break;
            k0 = k1;
            //cout << k0 << " ";
            if (k0 <= lb || k0 >= ub) {
                failed = true;
                break;
            }
            Eigen::MatrixXd der = asv(k0);
            *d1t = der(0, 0);
            *d2t = der(0, 1);
        }
        if (failed)
            std::cout << "not a local minimum" << std::endl;
        else {
            std::cout << "converged to " << k0 << " in " << iter << " iterations." << std::endl;
            loc_min.push_back(k0);
        }
    }

    toc = chrono::high_resolution_clock::now();

    std::cout << "Local minima of the smallest singular value:" << std::endl;
    for (std::vector<double>::const_iterator it = loc_min.begin(); it != loc_min.end(); ++it) {
        std::cout << *it << std::endl;
    }

    std::cout << "Total time: " << duration_cast<seconds>(toc - tic).count() << " sec" << std::endl;

    return 0;
}

