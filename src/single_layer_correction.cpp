/**
 * \file single_layer.cpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Single Layer BIO, using the transformations given in
 *        \f$\ref{ss:quadapprox}\f$ in the Lecture Notes for Advanced Numerical
 *        Methods for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "single_layer_correction.hpp"
#include "integral_gauss.hpp"

namespace parametricbem2d {
namespace single_layer_correction {
Eigen::MatrixXd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                  const AbstractParametrizedCurve &pi_p,
                                  const AbstractBEMSpace &space,
                                  const QuadRule &GaussQR,
                                  const double k) {
  double tol = std::numeric_limits<double>::epsilon();

  if (&pi == &pi_p) // Same Panels case
    return ComputeIntegralCoinciding(pi, pi_p, space, GaussQR,k);

  else if (fabs((pi(1) - pi_p(-1)).norm()) < tol ||
           fabs((pi(-1) - pi_p(1)).norm()) < tol) // Adjacent Panels case
    return ComputeIntegralAdjacent(pi, pi_p, space, GaussQR,k);

  else // Disjoint panels case
    return ComputeIntegralGeneral(pi, pi_p, space, GaussQR,k);
}

Eigen::MatrixXd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                          const AbstractParametrizedCurve &pi_p,
                                          const AbstractBEMSpace &space,
                                          const QuadRule &GaussQR,
                                          const double k) {
  unsigned N = GaussQR.n; // Quadrature order for the GaussQR object. Same order
                          // to be used for log weighted quadrature
  int Q = space.getQ(); // No. of Reference Shape Functions in trial/test space
  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // Lambda expression for functions F and G in \f$\eqref{eq:Vidp}\f$
      auto F = [&](double t) { // Function associated with panel pi_p
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s) { // Function associated with panel pi
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      // Lambda expression for the 1st integrand in \f$\eqref{eq:Isplit}\f$
      auto integrand1 = [&](double s, double t) {
        double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        double s_st;
        if (fabs(s - t) >
            sqrt_epsilon) // Away from singularity for stable evaluation as
                          // mentioned in \f$\eqref{eq:Stab}\f$
          // Simply evaluating the expression
          s_st = (pi(s) - pi_p(t)).squaredNorm() / (s - t) / (s - t);
        else // Near singularity
          // Using analytic limit for s - > t given in \f$\eqref{eq:Sdef}\f$
          s_st = pi.Derivative(0.5 * (t + s)).squaredNorm();
        return 0.5 * log(k*s_st) * F(t) * G(s);
      };

      double i1 = 0., i2 = 0.; // The two integrals in \f$\eqref{eq:Isplit}\f$

      // Tensor product quadrature for double 1st integral in
      // \f$\eqref{eq:Isplit}\f$
      for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
          i1 += GaussQR.w(i) * GaussQR.w(j) *
                integrand1(GaussQR.x(i), GaussQR.x(j));
        }
      }

      // Lambda expression for inner integrand in transformed coordinates in
      // \f$\eqref{eq:I21}\f$
      auto integrand2 = [&](double w, double z) {
        return F(0.5 * (w - z)/k) * G(0.5 * (w + z)/k) +
               F(0.5 * (w + z)/k) * G(0.5 * (w - z)/k);
      };

      // Getting log weighted quadrature nodes and weights
      QuadRule logweightQR = getLogWeightQR(k*2., N);

      // Double loop for 2nd double integral \f$\eqref{eq:I21}\f$
      for (unsigned int i = 0; i < N; ++i) {
        // Outer integral evaluated with Log weighted quadrature
        double z = logweightQR.x(i);
        double inner = 0.;
        // Evaluating the inner integral for fixed z
        for (unsigned int j = 0; j < N; ++j) {
          // Scaling Gauss Legendre quadrature nodes to the integral limits
          double w = GaussQR.x(j) * (2 - z);
          inner += GaussQR.w(j) * integrand2(w, z);
        }
        // Multiplying the integral with appropriate constants for
        // transformation to w from Gauss Legendre nodes
        inner *= (2 - z);
        i2 += logweightQR.w(i) * inner;
      }
      // Filling the matrix entry
      interaction_matrix(i, j) = -1. / (2. * M_PI) * (i1 + 0.5 * k*k * i2);
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                        const AbstractParametrizedCurve &pi_p,
                                        const AbstractBEMSpace &space,
                                        const QuadRule &GaussQR,
                                        const double k) {
  unsigned N = GaussQR.n; // Quadrature order for the GaussQR object. Same order
                          // to be used for log weighted quadrature
  int Q = space.getQ(); // No. of Reference Shape Functions in trial/test space
  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // Panel lengths for local arclength parametrization in
      // \f$\eqref{eq:ap}\f$. Actual values are not required so a length of 1 is
      // used for both the panels
      double length_pi = 1.;   // Length for panel pi
      double length_pi_p = 1.; // Length for panel pi_p

      // when transforming the parametrizations from [-1,1]->\Pi to local
      // arclength parametrizations [0,|\Pi|] -> \Pi, swap is used to ensure
      // that the common point between the panels corresponds to the parameter 0
      // in both arclength parametrizations
      bool swap = (pi(1) != pi_p(-1));

      // Lambda expressions for the functions F,G and D(r,phi) in
      // \f$\eqref{eq:Isplitapn}\f$
      auto F = [&](double t_pr) { // Function associated with panel pi_p
        // Transforming the local arclength parameter to standard parameter
        // range [-1,1] using swap
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s_pr) { // Function associated with panel pi
        // Transforming the local arclength parameter to standard parameter
        // range [-1,1] using swap
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      auto D_r_phi = [&](double r, double phi) { // \f$\eqref{eq:Ddef}\f$
        double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        // Transforming to local arclength parameter range
        double s_pr = r * cos(phi);
        // Transforming to standard parameter range [-1,1] using swap
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        // Transforming to local arclength parameter range
        double t_pr = r * sin(phi);
        // Transforming to standard parameter range [-1,1] using swap
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        if (r > sqrt_epsilon) // Away from singularity, simply use the formula
          return k*(pi(s) - pi_p(t)).squaredNorm() / r / r;
        else // Near singularity, use analytically evaluated limit at r -> 0 for
             // stable evaluation \f$\eqref{eq:Dstab}\f$
          return k*(1 - sin(2 * phi) * pi.Derivative(s).dot(pi_p.Derivative(t)));
      };

      // The two integrals in \f$\eqref{eq:Isplitapn}\f$ have to be further
      // split into two parts part 1 is where phi goes from 0 to alpha part 2 is
      // where phi goes from alpha to pi/2
      double alpha = atan(length_pi_p / length_pi); // the split point

      // i_IJ -> Integral I, part J
      double i11 = 0., i21 = 0., i12 = 0., i22 = 0.;
      // part 1 (phi from 0 to alpha)
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature node into phi
        double phi = alpha / 2 * (1 + GaussQR.x(i));
        // Computing inner integral with fixed phi
        // Inner integral for double integral 1, evaluated with Gauss Legendre
        // quadrature
        double inner1 = 0.;
        // Inner integral for double integral 2, evaluated with Log weighted
        // Gauss quadrature
        double inner2 = 0.;
        // Upper limit for inner 'r' integral
        double rmax = length_pi / cos(phi);
        // Evaluating the inner 'r' integral
        for (unsigned int j = 0; j < N; ++j) {
          // Getting Quadrature weights and nodes for Log weighted Gauss
          // quadrature
          QuadRule logweightQR = getLogWeightQR(rmax, N);
          // Evaluating inner2 using Log weighted Gauss quadrature
          double r = logweightQR.x(j);
          inner2 += logweightQR.w(j) * r * F(r * sin(phi)) * G(r * cos(phi));

          // Evaluating inner1 using Gauss Legendre quadrature
          r = rmax / 2 * (1 + GaussQR.x(j));
          inner1 += GaussQR.w(j) * r * log(D_r_phi(r, phi)) * F(r * sin(phi)) *
                    G(r * cos(phi));
        }
        // Multiplying the integral with appropriate constants for
        // transformation to r from Gauss Legendre nodes
        inner1 *= rmax / 2;
        // Multiplying the integrals with appropriate constants for
        // transformation to phi from Gauss Legendre nodes
        i11 += GaussQR.w(i) * inner1 * alpha / 2;
        i21 += GaussQR.w(i) * inner2 * alpha / 2;
      }

      // part 2 (phi from alpha to pi/2)
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature node into phi (alpha to pi/2)
        double phi =
            GaussQR.x(i) * (M_PI / 2. - alpha) / 2. + (M_PI / 2. + alpha) / 2.;
        // Computing inner integral with fixed phi
        // Inner integral for double integral 1, evaluated with Gauss Legendre
        // quadrature
        double inner1 = 0.;
        // Inner integral for double integral 2, evaluated with Log weighted
        // Gauss quadrature
        double inner2 = 0.;
        // Upper limit for inner 'r' integral
        double rmax = length_pi_p / sin(phi);
        // Evaluating the inner 'r' integral
        for (unsigned int j = 0; j < N; ++j) {
          // Getting Quadrature weights and nodes for Log weighted Gauss
          // quadrature
          QuadRule logweightQR = getLogWeightQR(rmax, N);
          // Evaluating inner2 using Log weighted Gauss quadrature
          double r = logweightQR.x(j);
          inner2 += logweightQR.w(j) * r * F(r * sin(phi)) * G(r * cos(phi));

          // Evaluating inner1 using Gauss Legendre quadrature
          r = rmax / 2 * (1 + GaussQR.x(j));
          inner1 += GaussQR.w(j) * r * log(D_r_phi(r, phi)) * F(r * sin(phi)) *
                    G(r * cos(phi));
        }
        // Multiplying the integral with appropriate constants for
        // transformation to r from Gauss Legendre quadrature nodes
        inner1 *= rmax / 2;
        // Multiplying the integrals with appropriate constants for
        // transformation to phi from Gauss Legendre quadrature nodes
        i12 += GaussQR.w(i) * inner1 * (M_PI / 2. - alpha) / 2.;
        i22 += GaussQR.w(i) * inner2 * (M_PI / 2. - alpha) / 2.;
      }
      // Summing up the parts to get the final integral
      double integral = 0.5 * (i11 + i12) + (i21 + i22);
      // Multiplying the integral with appropriate constants for transformation
      // to local arclength variables
      integral *= 4 / length_pi / length_pi_p;
      // Filling up the matrix entry
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                       const AbstractParametrizedCurve &pi_p,
                                       const AbstractBEMSpace &space,
                                       const QuadRule &GaussQR,
                                       const double k) {
  unsigned N = GaussQR.n; // Quadrature order for the GaussQR object.
  // Calculating the quadrature order for stable evaluation of integrands for
  // disjoint panels as mentioned in \f$\ref{par:distpan}\f$
  unsigned n0 = 30; // Order for admissible cases
  // Admissibility criteria
  double eta = 0.5;
  // Calculating the quadrature order
  unsigned n = n0 * std::max(1., 1. + log(rho(pi, pi_p) / eta));
  // std::cout << "ComputeIntegralGeneral used with order " << N << std::endl;
  // No. of Reference Shape Functions in trial/test space
  int Q = space.getQ();
  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // Lambda expression for functions F and G in \f$\eqref{eq:titg}\f$ for
      // Single Layer BIO
      auto F = [&](double t) { // Function associated with panel pi_p
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s) { // Function associated with panel pi
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      double integral = 0.;

      // Tensor product quadrature rule
      for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
          double s = GaussQR.x(i);
          double t = GaussQR.x(j);
          integral += GaussQR.w(i) * GaussQR.w(j) *
                      log(k*(pi(s) - pi_p(t)).norm()) * F(t) * G(s);
        }
      }
      // Filling up the matrix entry
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd GalerkinMatrix(const ParametrizedMesh mesh,
                               const AbstractBEMSpace &space,
                               const unsigned int &N,
                               const double k) {
  // Getting the number of panels in the mesh
  unsigned int numpanels = mesh.getNumPanels();
  // Getting dimensions of trial/test space
  unsigned int dims = space.getSpaceDim(numpanels);
  // Getting the panels from the mesh
  PanelVector panels = mesh.getPanels();
  // Getting the number of local shape functions in the trial/test space
  unsigned int Q = space.getQ();
  // Initializing the Galerkin matrix with zeros
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(dims, dims);
  // Panel oriented assembly \f$\ref{pc:ass}\f$
  QuadRule GaussQR = getGaussQR(N,-1.,1.);
  for (unsigned int i = 0; i < numpanels; ++i) {
    for (unsigned int j = 0; j < numpanels; ++j) {
      // Getting the interaction matrix for the pair of panels i and j
      Eigen::MatrixXd interaction_matrix =
          InteractionMatrix(*panels[i], *panels[j], space, GaussQR, k);
      // Local to global mapping of the elements in interaction matrix
      for (unsigned int I = 0; I < Q; ++I) {
        for (unsigned int J = 0; J < Q; ++J) {
          int II = space.LocGlobMap(I + 1, i + 1, numpanels) - 1;
          int JJ = space.LocGlobMap(J + 1, j + 1, numpanels) - 1;
          // Filling the Galerkin matrix entries
          output(II, JJ) += interaction_matrix(I, J);
        }
      }
    }
  }
  return output;
}

double Potential(const Eigen::Vector2d &x, const Eigen::VectorXd &coeffs,
                 const ParametrizedMesh &mesh, const AbstractBEMSpace &space,
                 const unsigned int &N) {
  // Getting the number of panels in the mesh
  unsigned int numpanels = mesh.getNumPanels();
  // Getting dimensions of BEM space
  unsigned int dims = space.getSpaceDim(numpanels);
  // asserting that the space dimension matches with coefficients
  assert(coeffs.rows() == dims);
  // Getting the panels from the mesh
  PanelVector panels = mesh.getPanels();
  // Getting the number of local shape functions in the BEM space
  unsigned int Q = space.getQ();
  // Getting general Gauss Quadrature rule
  QuadRule GaussQR = getGaussQR(N,-1.,1.);
  // Initializing the single layer potential vector for the bases, with zeros
  Eigen::VectorXd potentials = Eigen::VectorXd::Zero(dims);
  // Looping over all the panels
  for (unsigned panel = 0; panel < numpanels; ++panel) {
    // Looping over all reference shape functions
    for (unsigned i = 0; i < Q; ++i) {
      // The local integrand for a panel and a reference shape function
      auto integrand = [&](double t) {
        // Single Layer Potential
        return -1. / 2. / M_PI *
               log((x - panels[panel]->operator()(t)).norm()) *
               space.evaluateShapeFunction(i, t) *
               panels[panel]->Derivative(t).norm();
      };
      // Local to global mapping
      double local = ComputeIntegral(integrand, -1., 1., GaussQR);
      unsigned ii = space.LocGlobMap(i + 1, panel + 1, numpanels) - 1;
      // Filling the potentials vector
      potentials(ii) += local;
    }
  }
  return coeffs.dot(potentials);
}

} // namespace single_layer
} // namespace parametricbem2d
