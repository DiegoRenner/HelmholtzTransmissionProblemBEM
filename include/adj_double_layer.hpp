/**
 * \file adj_double_layer.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Adjoint Double Layer BIO, internally using the Double Layer
 *        implementation
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef ADJDOUBLELAYERHPP
#define ADJDOUBLELAYERHPP

#include <Eigen/Dense>
#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "parametrized_mesh.hpp"

namespace parametricbem2d {
/**
 * This namespace contains all the functions for evaluating the Adjoint Double
 * Layer Galerkin Matrix using quadrature and panel oriented assembly
 */
namespace adj_double_layer {

/**
 * This function is used to evaluate the full Galerkin matrix based on the
 * Bilinear form for Adjoint Double Layer BIO. It uses the trial and test spaces
 * and the parametrized mesh object, specified in the inputs to the function.
 * It evaluates the matrix by internally using the Double Layer Galerkin Matrix
 * evaluation and transposing it to give the final result.
 *
 * @param mesh ParametrizedMesh object containing all the panels in the form
 *             of small parametrized curves
 * @param trial_space The trial space for evaluating the matrix.
 * @param test_space The test space for evaluating the matrix.
 * @param N The order for gauss/log-weighted quadrature.
 * @return An Eigen::MatrixXd type Galerkin Matrix for the given mesh and space
 */
Eigen::MatrixXd GalerkinMatrix(const ParametrizedMesh mesh,
                               const AbstractBEMSpace &trial_space,
                               const AbstractBEMSpace &test_space,
                               const unsigned int &N);

} // namespace adj_double_layer
namespace adj_double_layer_helmholtz {
    Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                   const AbstractBEMSpace &trial_space,
                                   const AbstractBEMSpace &test_space,
                                   const unsigned int &N,
                                   const double c_o,
                                   const double c_i);


}
} // namespace parametricbem2d

#endif // ADJDOUBLELAYERHPP
