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

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "parametrized_mesh.hpp"

namespace parametricbem2d {
/**
 * This namespace contains all the functions for evaluating the Adjoint Double
 * Layer Galerkin Matrix using quadrature and panel oriented assembly
 */
namespace adj_double_layer_helmholtz {
    Eigen::MatrixXcd GalerkinMatrix(const ParametrizedMesh mesh,
                                   const AbstractBEMSpace &trial_space,
                                   const AbstractBEMSpace &test_space,
                                   const unsigned int &N,
                                   const double k);


} // namespace adj_double_layer_helmholtz
} // namespace parametricbem2d

#endif // ADJDOUBLELAYERHPP
