/**
 * \file logweight_quadrature.cpp
 * \brief This file defines the function to evaluate the quadrature rule
 *        for a log-weighted integral of the form :
 *        \f$ \int_{0}^{a} \log(t) f(t) dt \f$
 *        using the Gauss Laguerre Rule and a simple integral transformation.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "logweight_quadrature.hpp"

#include <Eigen/Dense>
#include "genLaguerreRule.hpp"

static const int KIND = 5; // To get Gauss-Laguerre quadrature rule

QuadRule getLogWeightQR(double a, int n){
  // get Gauss-Laguerre quadrature points and weights
  double *w = new double[n];
  double *x = new double[n];
  cgqf(n, KIND, 0, 0, 0, 1, x, w);
  // create new Quadrature rule
  QuadRule logWeightQR;
  logWeightQR.dim = 1;
  logWeightQR.n = n;
    std::cout << "test3" << std::endl;
    // create matrix and vector for its points and weights
  Eigen::MatrixXd points(n,1);
  Eigen::VectorXd weights(n);
  // Generate them applying the required transformation
  for(int i=0; i<n; i++){
    // quadrature point becomes: $q = a e^{-x}$
    points(i,0) = a*exp(-x[i]);
    // weight becomes $w = a(log(a)-x)w$
    weights(i) = a*(log(a)-x[i])*w[i];
  }
  logWeightQR.x = points;
  logWeightQR.w = weights;
  return logWeightQR;
}
