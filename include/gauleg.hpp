/**
 * \file gauleg.hpp
 * \brief This file defines the function to evaluate points and weights for
 *        Gauss Legendre quadrature rule.
 *
 * This File is a part of the 2D-Parametric BEM package
 */
#ifndef GAULEGHPP
#define GAULEGHPP

#include "logweight_quadrature.hpp"

/* Compute Gaussian quadrature nodes and weights for n nodes over
 * interval [a,b] \param[in] a,b Interval [a,b] endpoints \param[in] n Number of
 * quadrature points \param[out] xq,wq Pair of Gauss-Legendre quadrature points
 * and weights
 *
 * @param a Lower end of the domain
 * @param b Upper end of the domain
 * @param n Order for the quadrature rule
 * @param eps Tolerane level for the quadrature rule
 * @return std pair object containing weights and nodes in Eigen:VectorXd
 */
inline std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
gauleg(double a, double b, int n, double eps = 1.e-13) {

  /* Initialize weights and quadrature points */
  /*Eigen::VectorXd x(n);
  Eigen::VectorXd w(n);
  Eigen::VectorXd v(n*n);

  for(int i=0; i<n; i++){
    x[i] = 0.0;
  };

  for(int i=0; i<n-1; i++){
    w[i] = (i+1.0)/sqrt((2.0*i+1.0)*(2.0*i+3.0));
    //w[i] = (i+1.0)/sqrt((4.0*(i+1.)*(i+1.)-1.0));
  }
  int info;
  Eigen::MatrixXd  A(n,n);
  A.diagonal(-1)=w.segment(0,n-1);
  A.diagonal(1)=w.segment(0,n-1);
  A.diagonal(0)=x;

  Eigen::EigenSolver<Eigen::MatrixXd> es;
  es.compute(A,true);

  x = es.eigenvalues().real();
  A = es.eigenvectors().real();
  v = Eigen::Map<Eigen::VectorXd>(A.data(), n*n);
  // Compute weights
  for(int i=0; i<n; i++){
    x[i] = x[i];
    w[i] = 2*v[i*n]*v[i*n];
  }*/
   if (a > b)
      throw std::domain_error("Domain end points not ordered!");

    int i, j, m;
    double xmid, xlen, p1, p2, p3, dp1, z, z1, wqi;
    Eigen::RowVectorXd xq(n), wq(n);

    m = (n + 1) / 2;
    xmid = 0.5 * (a + b);
    xlen = 0.5 * (b - a);

    // get roots
    for (i = 0; i < m; i++) {

      // i-th root guess
      z = std::cos(M_PI * (i + 1 - 0.25) / (n + 0.5));

      // get i-th root
      do {
        p1 = 1.;
        p2 = 0.;
        for (j = 1; j <= n; j++) {
          p3 = p2;
          p2 = p1;
          p1 = ((2. * j - 1.) * z * p2 - (j - 1.) * p3) / j;
        }
        dp1 = n * (z * p1 - p2) / (z * z - 1.0);
        z1 = z;
        z = z1 - p1 / dp1;
      } while (std::abs(z - z1) > eps);

      // set nodes
      xq(i) = xmid - xlen * z;
      xq(n - 1 - i) = xmid + xlen * z;

      // set weights
      wqi = 2. * xlen / ((1. - z * z) * dp1 * dp1);
      wq(i) = wqi;
      wq(n - 1 - i) = wqi;
    }
  return std::make_pair(xq, wq);
}

inline std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
        cgauleg_redux(double a, double b, int n){
    int N = n*(n+1)/2-1;
    double sigma = (sqrt(2.)-1.)*(sqrt(2.)-1.);
    double xl = sigma;
    double xr = 1.;
    int ii = N-1;
    int qq = n;
    Eigen::RowVectorXd weights_gauleg, points_gauleg;
    Eigen::RowVectorXd weights(N), points(N);
    for (int i=1; i<n; i++){
        std::tie(points_gauleg, weights_gauleg) =
                gauleg(0, 1, qq, std::numeric_limits<double>::epsilon());
        for (int j=qq-1; j>=0; j--){
            weights[ii] = (xr-xl)*weights_gauleg[j];
            points[ii] = (xr-xl)*points_gauleg[j]+xl;
            ii--;
        };
        qq--;
        xr=xl;
        xl=xl*sigma;
    };
    return std::make_pair(points, weights);
}

/**
 * This function is evaluates a standard Gaussian Quadrature rule for the domain
 * [-1,1] for the given order. The quadrature rule is returned in the form of a
 * QuadRule object
 *
 * @param N Order for Gaussian Quadrature
 * @return QuadRule object containing the quadrature rule
 */
inline QuadRule getGaussQR(unsigned N) {
  // Getting standard Gauss Legendre Quadrature weights and nodes
  Eigen::RowVectorXd weights, points;
  std::tie(points, weights) =
      gauleg(0, 1, N, std::numeric_limits<double>::epsilon());
  QuadRule gauss;
  gauss.dim = 1;
  gauss.n = N;
  gauss.x = points;
  gauss.w = weights;
  return gauss;
}

inline QuadRule getCGaussQR(unsigned N) {
    // Getting standard Gauss Legendre Quadrature weights and nodes
    Eigen::RowVectorXd weights, points;
    std::tie(points, weights) =
            cgauleg_redux(0, 1, N);
    QuadRule gauss;
    gauss.dim = 1;
    gauss.n = N*(N+1)/2-1;
    gauss.x = points;
    gauss.w = weights;
    return gauss;
}
#endif
