/**
 * \file dirichlet.hpp
 * \brief This file defines lowest order indirect/direct BVP solvers to solve a
 * Dirichlet Boundary Value problem of the form given in \f$\eqref{eq:dirbvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef DIRICHLETHPP
#define DIRICHLETHPP

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "adj_double_layer.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "double_layer.hpp"
#include "hypersingular.hpp"
#include "integral_gauss.hpp"
#include "parametrized_mesh.hpp"
#include "single_layer.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace parametricbem2d {
/**
 * This function is used to evaluate Mass Matrices appearing in the variational
 * formulations of some boundary integral equations. The Mass matrix is of the
 * form \f$ M_{ij} = \int_{\Gamma} b^{i}(x) \beta^{j}(x) dS(x) \f$
 * Where \f$ b^{i} \f$ are the reference shape functions of test space and
 * \f$ \beta^{j} \f$ are the reference shape functions of the trial space.
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param x_space Space representing the rows / test space.
 * @param y_space Space representing the columns / trial space.
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::MatrixXd type representing the mass matrix
 */
inline Eigen::MatrixXd MassMatrix(const ParametrizedMesh &mesh,
                                  const AbstractBEMSpace &x_space,
                                  const AbstractBEMSpace &y_space,
                                  unsigned order) {
  // Getting the number of reference shape functions in both spaces
  unsigned qx = x_space.getQ();
  unsigned qy = y_space.getQ();
  // Getting the panels
  PanelVector panels = mesh.getPanels();
  unsigned numpanels = mesh.getNumPanels();
  unsigned split = numpanels/2;
  // Getting the space dimensions to fix matrix sizes
  unsigned int rows = x_space.getSpaceDim(numpanels);
  unsigned int cols = y_space.getSpaceDim(numpanels);
  // Setting the output matrix size
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows, cols);
  // Interaction matrix for all pairs of reference shape functions between the
  // two spaces
  Eigen::MatrixXd interaction(qx, qy);
  // Looping over all the panels
  for (unsigned panel = 0; panel < numpanels; ++panel) {
    for (unsigned i = 0; i < qx; ++i) {
      for (unsigned j = 0; j < qy; ++j) {
        // Transformed integrand for the mass matrix (in reference coordinates)
        std::function<double(double)> integrand = [&](double x) {
          return x_space.evaluateShapeFunction(i, x) *
                 y_space.evaluateShapeFunction(j, x) *
                 panels[panel]->Derivative(x).norm();
        };
        // Computing the interaction matrix using general Gauss quadrature
        interaction(i, j) = ComputeIntegral(integrand, -1, 1, order);
      }
    }
    // Local to global mapping of the elements in interaction
    for (unsigned int I = 0; I < qx; ++I) {
      for (unsigned int J = 0; J < qy; ++J) {
        // int II = x_space.LocGlobMap(I + 1, panel + 1, numpanels) - 1;
        // if(split == 0)
        // int JJ = y_space.LocGlobMap(J + 1, panel + 1, numpanels) - 1;
        int II = x_space.LocGlobMap2(I + 1, panel + 1, mesh) - 1;
        int JJ = y_space.LocGlobMap2(J + 1, panel + 1, mesh) - 1;
        // Filling the mass matrix entries
        output(II, JJ) += interaction(I, J);
      }
    }
  }
  return output;
}

/**
 * This namespace contains all the solvers for Dirichlet bvp of the form
 * \f$\eqref{eq:dirbvp}\f$. For different methods, the outputs mean something
 * else. Check the solver's documentation for further details.
 */


namespace dirichlet_bvp {
/**
 * This namespace contains the solver using the direct first kind method which
 * has the variational formulation as given in \f$\eqref{eq:aVdir}\f$. The
 * Solver uses the lowest order BEM spaces for computation.
 */
namespace direct_first_kind {
/**
 * This function is used to solve the Dirichlet boundary value problem given
 * in \f$\eqref{eq:dirbvp}\f$ using the variational formulation given in
 * \f$\eqref{eq:aVdir}\f$. The function outputs a vector of estimated Neumann
 * trace of the solution u.
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param g Dirichlet boundary condition in 2D using a function of the form
 *          double(double,double)
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd type representing the Neumann trace of the
 * solution u
 */
Eigen::VectorXd solve(const ParametrizedMesh &mesh,
                      std::function<double(double, double)> g, unsigned order) {
  // Same trial and test spaces
  DiscontinuousSpace<0> trial_space;
  DiscontinuousSpace<0> test_space;
  // Space used for interpolation of Dirichlet data
  ContinuousSpace<1> g_interpol_space;
  // Computing V matrix
  Eigen::MatrixXd V = single_layer::GalerkinMatrix(mesh, trial_space, order);
  // Computing K matrix
  Eigen::MatrixXd K =
      double_layer::GalerkinMatrix(mesh, g_interpol_space, test_space, order);
  // Computing mass matrix
  Eigen::MatrixXd M = MassMatrix(mesh, test_space, g_interpol_space, order);
  // Getting Dirichlet data at the vertices, interpolation by \f$S_{1}^{0}\f$
  Eigen::VectorXd g_N = g_interpol_space.Interpolate(g, mesh);
  // Build rhs for solving
  Eigen::VectorXd rhs = (0.5 * M + K) * g_N;
  // Solving for coefficients
  //Eigen::FullPivLU<Eigen::MatrixXd> dec(V);
  Eigen::HouseholderQR<Eigen::MatrixXd> dec(V);
  Eigen::VectorXd sol = dec.solve(rhs);
  return sol;
}
} // namespace direct_first_kind

/**
 * This namespace contains the solver using the direct second kind method which
 * has the variational formulation as given in \f$\eqref{eq:l2dv}\f$. The
 * Solver uses the lowest order BEM spaces for computation.
 */
namespace direct_second_kind {
/**
 * This function is used to solve the Dirichlet boundary value problem given
 * in \f$\eqref{eq:dirbvp}\f$ using the variational formulation given in
 * \f$\eqref{eq:l2dv}\f$ after lifting of the variation formulation given in
 * \f$\eqref{eq:bie2nbpvv}\f$ to the \f$L^{2}(\Gamma)\f$ space. The function
 * outputs a vector of estimated Neumann trace of the solution u.
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param g Dirichlet boundary condition in 2D using a function of the form
 *          double(double,double)
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd type representing the Neumann trace of the
 * solution u
 */
Eigen::VectorXd solve(const ParametrizedMesh &mesh,
                      std::function<double(double, double)> g, unsigned order) {
  // Same trial and test spaces
  ContinuousSpace<2> trial_space;
  ContinuousSpace<2> test_space;
  // Space used for interpolation of Dirichlet data
  ContinuousSpace<2> g_interpol_space;
  // Computing W matrix
  Eigen::MatrixXd W =
      hypersingular::GalerkinMatrix(mesh, g_interpol_space, order);
  // Computing K' matrix
  Eigen::MatrixXd Kp =
      adj_double_layer::GalerkinMatrix(mesh, trial_space, test_space, order);
  // Computing mass matrix
  Eigen::MatrixXd M = MassMatrix(mesh, test_space, trial_space, order);
  // Getting Dirichlet data
  Eigen::VectorXd g_N = g_interpol_space.Interpolate(g, mesh);
  // Build lhs for solving
  Eigen::MatrixXd lhs = (0.5 * M - Kp);
  // Build rhs for solving
  Eigen::VectorXd rhs = W * g_N;
  // Eigen::JacobiSVD<Eigen::MatrixXd> svd(lhs, ComputeThinU | ComputeThinV);
  // Eigen::VectorXd svals = svd.singularValues();
  // std::cout << "svals \n" << svals << std::endl;
  // std::cout << "Condition number: " << svals(0)/svals(M.rows()-1) <<
  // std::endl;
  Eigen::FullPivLU<Eigen::MatrixXd> dec(lhs);
  Eigen::VectorXd sol = dec.solve(rhs);
  return sol;
}
} // namespace direct_second_kind

/**
 * This namespace contains the solver using the indirect first kind method which
 * has the variational formulation as given in \f$\eqref{eq:iddirVv}\f$. The
 * Solver uses the lowest order BEM spaces for computation.
 */
namespace indirect_first_kind {
/**
 * This function is used to solve the Dirichlet boundary value problem given
 * in \f$\eqref{eq:dirbvp}\f$ using the variational formulation given in
 * \f$\eqref{eq:iddirVv}\f$. The function outputs a vector \f$\Phi\f$ such that
 * the solution u to the BVP can be constructed as \f$u =
 * \Psi^{\Delta}_{SL}(\Phi)\f$
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param g Dirichlet boundary condition in 2D using a function of the form
 *          double(double,double)
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd type representing \f$\Phi\f$ as described above
 */
Eigen::VectorXd solve(const ParametrizedMesh &mesh,
                      std::function<double(double, double)> g, unsigned order) {
  // Same trial and test spaces
  DiscontinuousSpace<0> trial_space;
  DiscontinuousSpace<0> test_space;
  // Space used for interpolation of Dirichlet data
  ContinuousSpace<1> g_interpol_space;
  // Computing V matrix
  Eigen::MatrixXd V = single_layer::GalerkinMatrix(mesh, trial_space, order);
  // Computing mass matrix
  Eigen::MatrixXd M = MassMatrix(mesh, test_space, g_interpol_space, order);
  // Getting Dirichlet data at the vertices, interpolation by \f$S_{1}^{0}\f$
  Eigen::VectorXd g_N = g_interpol_space.Interpolate(g, mesh);
  // Build rhs for solving
  Eigen::VectorXd rhs = M * g_N;
  // Solving for coefficients
  Eigen::FullPivLU<Eigen::MatrixXd> dec(V);
  Eigen::VectorXd sol = dec.solve(rhs);
  return sol;
}
} // namespace indirect_first_kind

/**
 * This namespace contains the solver using the indirect second kind method. The
 * Solver uses the lowest order BEM spaces for computation.
 */
namespace indirect_second_kind {
/**
 * This function is used to solve the Dirichlet boundary value problem given
 * in \f$\eqref{eq:dirbvp}\f$ using the variational formulation of indirect
 * second kind BIE. The function outputs a vector \f$\Phi\f$ such that the
 * solution u to the BVP can be constructed as \f$u =
 * \Psi^{\Delta}_{DL}(\Phi)\f$
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param g Dirichlet boundary condition in 2D using a function of the form
 *          double(double,double)
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd type representing \f$\Phi\f$ as described above
 */
Eigen::VectorXd solve(const ParametrizedMesh &mesh,
                      std::function<double(double, double)> g, unsigned order) {
  // Same trial and test spaces
  DiscontinuousSpace<0> trial_space;
  DiscontinuousSpace<0> test_space;
  // Space used for interpolation of Dirichlet data
  ContinuousSpace<1> g_interpol_space;
  // Computing K matrix
  Eigen::MatrixXd K =
      double_layer::GalerkinMatrix(mesh, trial_space, trial_space, order);
  // Computing mass matrix for lhs
  Eigen::MatrixXd Ml = MassMatrix(mesh, test_space, trial_space, order);
  // Computing mass matrix for rhs
  Eigen::MatrixXd Mr = MassMatrix(mesh, test_space, g_interpol_space, order);
  // Getting Dirichlet data at the vertices, interpolation by \f$S_{1}^{0}\f$
  Eigen::VectorXd g_N = g_interpol_space.Interpolate(g, mesh);
  // Build lhs for solving
  Eigen::MatrixXd lhs = (-0.5 * Ml + K);
  // Build rhs for solving
  Eigen::VectorXd rhs = Mr * g_N;
  // Solving for coefficients
  Eigen::FullPivLU<Eigen::MatrixXd> dec(lhs);
  Eigen::VectorXd sol = dec.solve(rhs);
  return sol;
}
} // namespace indirect_second_kind
} // namespace dirichlet_bvp
} // namespace parametricbem2d

#endif // DIRICHLETHPP
