/**
 * \file roots_brent_square_arnoldi.cpp
 * \brief This target builds a script that computes minimas in 
 * the smallest singular value of the
 * Galerkin BEM approximated solutions operator for the sedond-kind direct 
 * BIEs of the Helmholtz
 * transmission problem using the Van Wijngaarden-Dekker-Brent method.
 * The scatterer is set to be a square.
 * The results are written to disk.
 * The script can be run as follows:
 *
 * <tt>
 *  /path/to/roots_brent_circle \<half side length of square\> \<refraction inside\>
 *     \<refraction outside\> \<initial wavenumber\> \<\#grid points for root search\>
 *     \<\#panels\> \<order of quadrature rule\> \<outputfile\>.
 * </tt>
 *
 * The resulting file will contain the left boundary of the 
 * interval used to compute the root in the first column. 
 * Then in the next three columns will be the point, the 
 * function value and the derivative at which the root was found.
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
inline double rsv(const Eigen::MatrixXcd &T, const Eigen::MatrixXcd &W, int q) {
    int nr = W.rows(), nc = W.cols();
    Eigen::MatrixXcd U, V, Q, B, C, thinQ = Eigen::MatrixXcd::Identity(nr, nc);
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu_decomp(T);
    V = lu_decomp.solve(W);
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr(V);
    Q = qr.householderQ() * thinQ;
    for (int i = 0; i < q; ++i) {
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr1(lu_decomp.adjoint().solve(Q));
        Q = qr1.householderQ() * thinQ;
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr2(lu_decomp.solve(Q));
        Q = qr2.householderQ() * thinQ;
    }
    B = lu_decomp.solve(Q);
    C = Q.transpose().conjugate() * B;
    Eigen::BDCSVD<Eigen::MatrixXcd> svd(C * Q.transpose().conjugate());
    return 1.0 / svd.singularValues()(0);
}

// Brent's algorithm routine
double local_min_rc (double &a, double &b, int &status, double value, double epsilon);

// Entry point
int main(int argc, char** argv) {

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
    double acc = atof(argv[8]);

    // generate output filename with set parameters
    std::string base_name = "../data/file_roots_brent_square_arnoldi_";
    std::string suffix = ".dat";
    std::string divider = "_";
    std::string file_minimas = base_name.append(argv[2]).append(divider).append(argv[5])
                                       .append(divider).append(argv[8]) + suffix;
    // clear existing file
    std::ofstream file_out;
    file_out.open(file_minimas, std::ofstream::out | std::ofstream::trunc);
    file_out.close();

    int nc = 2, q = 2;
    int nr = 2 * numpanels;
    Eigen::MatrixXcd W = standard_normal_random_matrix(nr, nc);

    double k_max = 4.0, k_step = (k_max - k_min) / n_points_k;

    std::cout << "Bracketing local minima...    [--------------------]";

    auto tic = chrono::high_resolution_clock::now();

    unsigned ksize = 1 + (int)round((k_max - k_min) / k_step), finished = 0;
    std::vector<double> kv(ksize), res(ksize), loc_min, bracket_left, bracket_right;
    std::generate(kv.begin(), kv.end(), [&]() { static double k = k_min - k_step; return k += k_step; });

    std::transform(std::execution::par_unseq, kv.cbegin(), kv.cend(), res.begin(), [&](double k) {
        Eigen::MatrixXcd T = gen_sol_op(mesh, order, k, c_o, c_i);
        double sv = rsv(T, W, q);
        ++finished;
        size_t progress = (100 * finished) / ksize, l1 = progress / 5, l2 = 20 - l1;
        std::cout << "\rBracketing local minima...    [" << std::string(l1, '=') << std::string(l2, '-') << "] (" << progress << "%)";
        std::flush(std::cout);
        return sv;
    });

    std::cout << "\rBracketing local minima...    Done.                          " << std::endl;

    bool lb = false;

    for (size_t i = 0; i < res.size() - 2; ++i) {
        double &c = res[i+1], L = c - res[i], R = res[i+2] - c;
        double k = k_min + i * k_step;
        if (!lb && L < 0. && R > 0.) { // local minimum
            bracket_left.push_back(k);
        }
        else if (L > 0. && R < 0.) { // local maximum
            if (lb) {
                bracket_right.push_back(k);
                lb = false;
            } else if (bracket_left.empty()) {
                bracket_left.push_back(k);
                lb = true;
            } else bracket_right.push_back(k);
        }
    }
    if (bracket_right.size() < bracket_left.size())
        bracket_right.push_back(k_max);
    assert (bracket_left.size() == bracket_right.size());

    auto toc = chrono::high_resolution_clock::now();

    size_t loc_min_count = bracket_left.size(); // number of local minima

    std::cout << "Found approximate locations of " << loc_min_count << " local minima." << std::endl;
    std::cout << "Elapsed time: " << duration_cast<seconds>(toc - tic).count() << " sec" << std::endl;

    // Routine for computing the smallest singular value by Arnoldi iterations
    auto asv = [&](double k) {
        Eigen::MatrixXcd T = gen_sol_op(mesh, order, k, c_o, c_i);
        return arnoldi::sv(T, 1, acc)(0);
    };

    for (size_t i = 0; i < loc_min_count; ++i) {
        std::cout << "Finding local minimum " << i + 1 << "...    ";
        std::flush(std::cout);
        double a = bracket_left[i], b = bracket_right[i], h = b - a, lp = a;
        int status = 0, iter = 0;
        double arg = local_min_rc(a, b, status, 0., epsilon_fin);
        while (status) {
            arg = local_min_rc(a, b, status, asv(arg), epsilon_fin);
            ++iter;
        }
        if (arg > lp + h * .95) {
            std::cout << "Not a local minimum" << std::endl;
        } else {
            std::cout << "Converged in " << iter << " iterations." << std::endl;
            loc_min.push_back(arg);
        }
    }

    toc = chrono::high_resolution_clock::now();

    std::cout << "Local minima of the smallest singular value:" << std::endl;
    for (std::vector<double>::const_iterator it = loc_min.begin(); it != loc_min.end(); ++it) {
        std::cout << *it << std::endl;
    }

    std::cout << "Total time: " << duration_cast<seconds>(toc - tic).count() << " sec" << std::endl;

    // Curves
#if 1
    std::cout << "rsv = [";
    for (std::vector<double>::const_iterator it = res.begin(); it != res.end(); ++it) {
        std::cout << *it << (it + 1 == res.end() ? "]" : " ");
    }
    std::cout << std::endl << "asv = [";
    for (double k = k_min; k <= k_max; k += k_step) {
        std::cout << asv(k) << (k == k_max ? "]" : " ");
    }
    std::cout << std::endl;
#endif

    // Timing
#if 0
    for (double k = k_min; k <= k_max; k += k_step) {
        Eigen::MatrixXcd T = gen_sol_op(mesh, order, k, c_o, c_i);
        tic = chrono::high_resolution_clock::now();
        arnoldi::sv(T, 1, acc)(0);
        toc = chrono::high_resolution_clock::now();
        auto dur1 = duration_cast<microseconds>(toc - tic);
        tic = chrono::high_resolution_clock::now();
        rsv(T, W, q);
        toc = chrono::high_resolution_clock::now();
        auto dur2 = duration_cast<microseconds>(toc - tic);
        std::cout << dur1.count() << " " << dur2.count() << "; " << std::endl;
    }
    std::cout << std::endl;
#endif

    return 0;
}

//****************************************************************************80

double local_min_rc ( double &a, double &b, int &status, double value, double epsilon )

//****************************************************************************80
//
//  Purpose:
//
//    LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
//
//  Discussion:
//
//    This routine seeks an approximation to the point where a function
//    F attains a minimum on the interval (A,B).
//
//    The method used is a combination of golden section search and
//    successive parabolic interpolation.  Convergence is never much
//    slower than that for a Fibonacci search.  If F has a continuous
//    second derivative which is positive at the minimum (which is not
//    at A or B), then convergence is superlinear, and usually of the
//    order of about 1.324...
//
//    The routine is a revised version of the Brent local minimization
//    algorithm, using reverse communication.
//
//    It is worth stating explicitly that this routine will NOT be
//    able to detect a minimizer that occurs at either initial endpoint
//    A or B.  If this is a concern to the user, then the user must
//    either ensure that the initial interval is larger, or to check
//    the function value at the returned minimizer against the values
//    at either endpoint.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters
//
//    Input/output, double &A, &B.  On input, the left and right
//    endpoints of the initial interval.  On output, the lower and upper
//    bounds for an interval containing the minimizer.  It is required
//    that A < B.
//
//    Input/output, int &STATUS, used to communicate between
//    the user and the routine.  The user only sets STATUS to zero on the first
//    call, to indicate that this is a startup call.  The routine returns STATUS
//    positive to request that the function be evaluated at ARG, or returns
//    STATUS as 0, to indicate that the iteration is complete and that
//    ARG is the estimated minimizer.
//
//    Input, double VALUE, the function value at ARG, as requested
//    by the routine on the previous call.
//
//    Output, double LOCAL_MIN_RC, the currently considered point.
//    On return with STATUS positive, the user is requested to evaluate the
//    function at this point, and return the value in VALUE.  On return with
//    STATUS zero, this is the routine's estimate for the function minimizer.
//
//  Local parameters:
//
//    C is the squared inverse of the golden ratio.
//
//    EPS is the square root of the relative machine precision.
//
{
  static double arg;
  static double c;
  static double d;
  static double e;
  static double eps;
  static double fu;
  static double fv;
  static double fw;
  static double fx;
  static double midpoint;
  static double p;
  static double q;
  static double r;
  static double tol;
  static double tol1;
  static double tol2;
  static double u;
  static double v;
  static double w;
  static double x;
//
//  STATUS (INPUT) = 0, startup.
//
  if ( status == 0 )
  {
    if ( b <= a )
    {
      cout << "\n";
      cout << "LOCAL_MIN_RC - Fatal error!\n";
      cout << "  A < B is required, but\n";
      cout << "  A = " << a << "\n";
      cout << "  B = " << b << "\n";
      status = -1;
      exit ( 1 );
    }
    c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

    eps = sqrt (epsilon);
    tol = epsilon;

    v = a + c * ( b - a );
    w = v;
    x = v;
    e = 0.0;

    status = 1;
    arg = x;

    return arg;
  }
//
//  STATUS (INPUT) = 1, return with initial function value of FX.
//
  else if ( status == 1 )
  {
    fx = value;
    fv = fx;
    fw = fx;
  }
//
//  STATUS (INPUT) = 2 or more, update the data.
//
  else if ( 2 <= status )
  {
    fu = value;

    if ( fu <= fx )
    {
      if ( x <= u )
      {
        a = x;
      }
      else
      {
        b = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x )
      {
        a = u;
      }
      else
      {
        b = u;
      }

      if ( fu <= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }
//
//  Take the next step.
//
  midpoint = 0.5 * ( a + b );
  tol1 = eps * std::abs ( x ) + tol / 3.0;
  tol2 = 2.0 * tol1;
//
//  If the stopping criterion is satisfied, we can exit.
//
  if ( fabs ( x - midpoint ) <= ( tol2 - 0.5 * ( b - a ) ) )
  {
    status = 0;
    return arg;
  }
//
//  Is golden-section necessary?
//
  if ( std::abs ( e ) <= tol1 )
  {
    if ( midpoint <= x )
    {
      e = a - x;
    }
    else
    {
      e = b - x;
    }
    d = c * e;
  }
//
//  Consider fitting a parabola.
//
  else
  {
    r = ( x - w ) * ( fx - fv );
    q = ( x - v ) * ( fx - fw );
    p = ( x - v ) * q - ( x - w ) * r;
    q = 2.0 * ( q - r );
    if ( 0.0 < q )
    {
      p = - p;
    }
    q = std::abs ( q );
    r = e;
    e = d;
//
//  Choose a golden-section step if the parabola is not advised.
//
    if (
      ( std::abs ( 0.5 * q * r ) <= std::abs ( p ) ) ||
      ( p <= q * ( a - x ) ) ||
      ( q * ( b - x ) <= p ) )
    {
      if ( midpoint <= x )
      {
        e = a - x;
      }
      else
      {
        e = b - x;
      }
      d = c * e;
    }
//
//  Choose a parabolic interpolation step.
//
    else
    {
      d = p / q;
      u = x + d;

      if ( ( u - a ) < tol2 )
      {
        d = tol1 * ( midpoint - x < 0. ? -1.0 : 1.0 );
      }

      if ( ( b - u ) < tol2 )
      {
        d = tol1 * ( midpoint - x < 0. ? -1.0 : 1.0 );
      }
    }
  }
//
//  F must not be evaluated too close to X.
//
  if ( tol1 <= std::abs ( d ) )
  {
    u = x + d;
  }
  if ( std::abs ( d ) < tol1 )
  {
    u = x + tol1 * ( d < 0. ? -1.0 : 1.0 );
  }
//
//  Request value of F(U).
//
  arg = u;
  status = status + 1;

  return arg;
}
