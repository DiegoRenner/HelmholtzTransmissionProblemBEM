/**
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef INTEGRALGAUSSHPP
#define INTEGRALGAUSSHPP

#include <limits>
#include <vector>

#include "gauleg.hpp"
#include "genLaguerreRule.hpp"
#include <Eigen/Dense>
namespace parametricbem2d {
/* This function computes an integral numerically using Gauss Legendre
 * Quadrature Rule of the given order
 *
 * @tparam T Template type for integrand. Should support evaluation.
 * @param integrand The integrand to be integrated
 * @param a Lower end of the integration domain
 * @param b Upper end of the integration domain
 * @param N Order for the quadrature rule
 * @return Integral value
 */
template <typename T>
double ComputeIntegral(T integrand, double a, double b, unsigned int N) {
  if (b < a) {
    throw std::domain_error("Domain end points not ordered!");
  }
  // Getting quadrature weights and points
  Eigen::RowVectorXd weights, points;
  std::tie(points, weights) =
      gauleg(a, b, N, std::numeric_limits<double>::epsilon());
  double integral = 0.;
  for (unsigned int i = 0; i < N; ++i)
    integral += weights(i) * integrand(points(i));
  return integral;
}

/* This function computes an integral numerically using Gauss Legendre
 * Quadrature Quadrature Rule of the given order
 *
 * @tparam T Template type for integrand. Should support evaluation.
 * @param integrand The integrand to be integrated
 * @param a Lower end of the integration domain
 * @param b Upper end of the integration domain
 * @param GaussQR Gaussian Quadrature in form of a QuadRule object
 * @return Integral value
 */
template <typename T>
double ComputeIntegral(T integrand, double a, double b,
                       const QuadRule &GaussQR) {
  if (b < a) {
    throw std::domain_error("Domain end points not ordered!");
  }
  // Getting quadrature weights and points
  double integral = 0.;
  unsigned N = GaussQR.n;
  double mean = 0.5 * (b + a);
  double diff = 0.5 * (b - a);
  for (unsigned int i = 0; i < N; ++i) {
    double x = GaussQR.x(i) * diff + mean;
    integral += GaussQR.w(i) * integrand(x);
  }
  return integral * diff;
}

/* This function computes a log weighted integral numerically using a Log
 * weighted Quadrature Rule of the given order, derived from Gauss Laguerre rule
 *
 * @tparam T Template type for integrand. Should support evaluation.
 * @param integrand The integrand to be integrated (without the log weight)
 * @param a Upper end of the integration domain
 * @param N Order for the quadrature rule
 * @return Integral value
 */
/*template <typename T>
double ComputeLogIntegral(T integrand, double a, unsigned int N) {
  QuadRule logweightQR = getLogWeightQR(a, N);
  double integral = 0.;
  for (unsigned int i = 0; i < N; ++i) {
    double x = logweightQR.x(i);
    integral += logweightQR.w(i) * integrand(x);
  }
  return integral;
}*/

/* This function computes a log weighted integral numerically using a Log
 * weighted Quadrature Rule of the given order, derived as Generalized Gauss
 * Quadrature with arbitrary weight function. This function is valid only for
 * integration in the following range: [0,1]
 *
 * @tparam T Template type for integrand. Should support evaluation.
 * @param integrand The integrand to be integrated (without the log weight)
 * @param N Order for the quadrature rule
 * @return Integral value
 */

/* This function computes a log weighted integral numerically using a Log
 * weighted Quadrature Rule of the given order. The integration domain [0.a] is
 * transformed to [0,1] and It is evaluated using Gauss Legendre rule and the
 * Log Weighted Quadrature rule which is derive as a Generalized Gauss
 * Quadrature rule for an arbitrary weight function.
 *
 * @tparam T Template type for integrand. Should support evaluation.
 * @param integrand The integrand to be integrated (without the log weight)
 * @param a Upper end of the integration domain
 * @param QR Gaussian Quadrature in form of a QuadRule object
 * @return Integral value
 */

/* This function computes a Gauss Lguerre integral using the respective
 * Quadrature Rule.
 *
 * @tparam T Template type for integrand. Should support evaluation.
 * @param integrand The integrand to be integrated (without the log weight)
 * @param a Upper end of the integration domain
 * @param n Order for the quadrature rule
 * @return Integral value
 */
template <typename T>
double ComputeLaguerreIntegral(T integrand, unsigned int n) {
  int KIND = 5; // To get Gauss-Laguerre quadrature rule
  double *w = new double[n];
  double *x = new double[n];
  cgqf(n, KIND, 0, 0, 0, 1, x, w);
  double integral = 0.;
  for (unsigned int i = 0; i < n; ++i) {
    integral += w[i] * integrand(x[i]);
  }
  delete[] w;
  delete[] x;
  return integral;
}

/* This function computes a double integral numerically using Gauss Legendre
 * Quadrature Quadrature Rule of the given order. The integral has to be of the
 * form \f$ \int_{x=a}^{b} \int_{y=ll(x)}^{ul(x)} f(x,y) dy dx \f$
 *
 * @tparam F Template type for integrand. Should support evaluation of the form
             f(x,y)
 * @param integrand The integrand to be integrated: f(x,y)
 * @param a Lower end of the integration domain for the variable x
 * @param b Upper end of the integration domain for the variable y
 * @tparam UL Template type for the upper limit function for the variable y in
              the double integral
 * @tparam LL Template type for the lower limit function for the variable y in
              the double integral
 * @param ll Function representing the lower limit of the y integral, as a
             function of x
 * @param ul Function representing the upper limit of the y integral, as a
             function of x
 * @param GaussQR Gaussian Quadrature in form of a QuadRule object
 * @return Integral value
 */
template <typename F, typename LL, typename UL>
double ComputeDoubleIntegral(const F &f, double a, double b, const LL &ll,
                             const UL &ul, unsigned order) {
  double integral = 0;
  // The function fx(x) obtained by evaluating the inner integral over y.
  auto fx = [&](double x) {
    // Lambda function to fix the x value in f(x,y) for integration over y
    auto temp = [&](double y) { return f(x, y); };
    // The limits for the inner y integral
    double A = ll(x); // lower limit
    double B = ul(x); // upper limit
    // Computing the inner integral over y
    return parametricbem2d::ComputeIntegral(temp, A, B, order);
  };
  // Integrating fx(x) to obtain the double integral value
  return parametricbem2d::ComputeIntegral(fx, a, b, order);
}

} // namespace parametricbem2d

#endif // INTEGRALGAUSSHPP
