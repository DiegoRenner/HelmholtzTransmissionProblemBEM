/**
 * \file gauleg.hpp
 * \brief This file defines the function to evaluate points and weights for
 *        Gauss Legendre quadrature rule.
 *
 * This File is a part of the HelmholtzTransmissionBEM.
 * It has been adapted from the 2D-Parametric BEM package.
 */

#ifndef GAULEGHPP
#define GAULEGHPP

/**
 * This Struct object is used to store a quadrature Rule.
 */
struct QuadRule {
    std::size_t dim;   // dimension of space
    std::size_t n;     // number of nodes/weights
    Eigen::MatrixXd x; // quadrature nodes (columns of a matrix with dim rows)
    Eigen::VectorXd w; // vector of quadrature weights
};

/**
 * Compute Gaussian quadrature nodes and weights for n nodes over
 * interval [a,b]. The nodes and weights are returned as a std pair.
 *
 * @param a Lower end of the domain
 * @param b Upper end of the domain
 * @param n Order for the quadrature rule
 * @param eps Tolerane level for the quadrature rule
 * @return std pair object containing weights and nodes in Eigen:VectorXd
 */
inline std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
gauleg(double a, double b, int n, double eps = 1.e-13) {

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

/**
 * Compute reduced composite Gaussian quadrature nodes and weights for n*(n+1)/2-1 nodes over
 * interval [0,1]. The nodes and weights are returned as a std pair.
 *
 * @param n Order for the quadrature rule
 * @param eps Tolerane level for the quadrature rule
 * @return std pair object containing weights and nodes in Eigen:VectorXd
 */
inline std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
        cgauleg_redux(int n, double eps){
    // total number of points
    int N = n*(n+1)/2-1;
    // set constant for geometric grading
    double sigma = (sqrt(2.)-1.)*(sqrt(2.)-1.);
    // initialize boundaries for first gauss-legendre rule
    double xl = sigma;
    double xr = 1.;
    // counters for mapping from gauss-legendre
    // to composite gauss-legendre
    int ii = N-1;
    int qq = n;
    // initializing points and weights
    Eigen::RowVectorXd weights_gauleg, points_gauleg;
    Eigen::RowVectorXd weights(N), points(N);
    for (int i=1; i<n; i++){
        // generate gauss-legendre rule for Interval [0,1]
        std::tie(points_gauleg, weights_gauleg) =
                gauleg(0, 1, qq, eps);
        for (int j=qq-1; j>=0; j--){
            // map to current boundaries and assign weights and points
            weights[ii] = (xr-xl)*weights_gauleg[j];
            points[ii] = (xr-xl)*points_gauleg[j]+xl;
            ii--;
        };
        qq--;
        // set new boundaries
        xr=xl;
        xl=xl*sigma;
    };
    return std::make_pair(points, weights);
}

/**
 * This function is evaluates a standard Gaussian quadrature rule for the domain
 * [a,b] for the given order. The quadrature rule is returned in the form of a
 * QuadRule object
 *
 * @param N Order for Gaussian quadrature
 * @param a left boundary of domain
 * @param b right boundary of domain
 * @return QuadRule object containing the quadrature rule
 */
inline QuadRule getGaussQR(unsigned N, double a, double b) {
  // Getting standard Gauss Legendre quadrature weights and nodes
  Eigen::RowVectorXd weights, points;
  std::tie(points, weights) =
      gauleg(a, b, N, std::numeric_limits<double>::epsilon());
  // assign to QuadRule struct
  QuadRule gauss;
  gauss.dim = 1;
  gauss.n = N;
  gauss.x = points;
  gauss.w = weights;
  return gauss;
}

/**
 * This function is evaluates a composite Gaussian quadrature rule for the domain
 * [0,1] for the given order. The quadrature rule is returned in the form of a
 * QuadRule object
 *
 * @param N Order for Gaussian quadrature
 * @return QuadRule object containing the quadrature rule
 */
inline QuadRule getCGaussQR(unsigned N) {
    // getting composite Gauss-Legendre quadrature weights and nodes
    Eigen::RowVectorXd weights, points;
    std::tie(points, weights) =
            cgauleg_redux(N, std::numeric_limits<double>::epsilon());
    // assign to QuadRule struct
    QuadRule gauss;
    gauss.dim = 1;
    gauss.n = N*(N+1)/2-1;
    gauss.x = points;
    gauss.w = weights;
    return gauss;
}
#endif
