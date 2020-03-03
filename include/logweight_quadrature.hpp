/**
 * \file logweight_quadrature.hpp
 * \brief This file declares the functions to evaluate the quadrature rule
 *        for a log-weighted integral of the form :
 *        \f$ \int_{0}^{a} \log(t) f(t) dt \f$
 *        using the Gauss Laguerre Rule and a simple integral transformation.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef LOGWEIGHTGQ
#define LOGWEIGHTGQ

#include <Eigen/Dense>
#include <utility>
#include <vector>

/**
 * This Struct object is used to store a Quadrature Rule.
 */
struct QuadRule {
  std::size_t dim;   // dimension of space
  std::size_t n;     // number of nodes/weights
  Eigen::MatrixXd x; // quadrature nodes (columns of a matrix with dim rows)
  Eigen::VectorXd w; // vector of quadrature weights
};

class test {
public:
  unsigned x;
  static std::vector<double> testvec;
  // public:
  static std::vector<double> returntestvec();
  unsigned wtfunction();
  static void init();
};

/**
 * This function evaluates the Log weighted quadrature rule and stores the
 * result in the struct QuadRule. The quadrature rule is computed using the
 * Gauss Laguerre rule by transforming the log-weighted integral to the
 * following form : \f$ a\int_{0}^{\infty} e^{-s} f(a e^{-s})(\log(a)-s) ds \f$
 *
 * @param a The upper limit for the integral
 * @param n Desired order for the quadrature rule
 * @return Quadrature rule in the form of QuadRule struct
 */
QuadRule getLogWeightQR(double a, int n);

/**
 * This function returns the log weighted quadrature rule derived from the
 * Mathematica script provided by * \f$\ref{}\f$ used to evaluate Generalized
 * Gaussian quadrature rules for arbitrary weight functions. The Quadrature rule
 * is stored and is available only at orders equal to powers of 2, until 256.
 * They are intended for the integration domain [0,1]
 *
 * @param N Desired order for the quadrature rule
 * @return std pair of vectors, containing points and weights of the quadrature
 * rule
 */
std::pair<std::vector<double>, std::vector<double>> getLogWeightQR(int N);

#endif
